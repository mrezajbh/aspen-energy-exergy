"""Human-readable tables and a machine-readable snapshot of one analysis."""

import csv
import io
import json
import os
from collections import defaultdict
from pathlib import Path
from tempfile import NamedTemporaryFile
from zipfile import ZIP_DEFLATED, ZipFile

from .quality import QUALITY_COLUMNS, build_quality


def _fmt(value, digits=3):
    return "unavailable" if value is None else f"{value:,.{digits}f}"


def summary(report):
    """Format the plant boundary and point to any incomplete or signed results."""
    plant, coverage = report["plant"], report["coverage"]
    lines = [
        "Aspen Plus energy and exergy analysis",
        f"Source: {report['source'].get('source_filename', 'normalized Aspen export')}",
        f"Settings: {report['settings'].get('name', 'user settings' if report['settings'] else 'as supplied')}",
        "Environment: 298.15 K, 101325 Pa. Steady state. Heat and work positive inward.",
        "",
        f"Material streams: {coverage['streams']}    Equipment blocks: {coverage['blocks']}",
        f"Stream total exergies: {coverage['total_exergy_available']}/{coverage['streams']}",
        f"Physical and chemical splits: {coverage['splits_available']}/{coverage['streams']}",
        "Split quality: " + ", ".join(f"{k} {v}" for k, v in sorted(coverage["split_quality"].items())),
        f"Equipment exergy balances: {coverage['block_exergy_balances_available']}/{coverage['blocks']}",
        "",
        "External plant boundary",
    ]
    for label, key, unit in [
        ("Material enthalpy in", "enthalpy_in_kw", "kW"),
        ("Material enthalpy out", "enthalpy_out_kw", "kW"),
        ("Heat in", "heat_input_kw", "kW"),
        ("Heat rejected", "heat_rejected_kw", "kW"),
        ("Work in", "work_input_kw", "kW"),
        ("Work out", "work_output_kw", "kW"),
        ("First-law residual", "energy_residual_kw", "kW"),
        ("Exergy in", "exergy_input_kw", "kW"),
        ("Exergy out", "exergy_output_kw", "kW"),
        ("Signed exergy destruction", "exergy_destruction_kw", "kW"),
        ("Entropy generation", "entropy_generation_kw_k", "kW/K"),
    ]:
        digits = 6 if unit == "kW/K" else 3
        lines.append(f"{label:<28} {_fmt(plant[key], digits):>16} {unit}")
    if plant["product_streams"]:
        lines.extend(["", "Selected material products"])
        for stream_id in plant["product_streams"]:
            stream = report["streams"][stream_id]
            lines.append(f"{stream_id:<16} {_fmt(stream.get('mass_kg_h')):>14} kg/h"
                         f"  {_fmt(stream['total_exergy_kw']):>14} kW total exergy")
        fraction = plant["product_exergy_fraction_of_input"]
        pct = None if fraction is None else 100.0 * fraction
        lines.append(f"Selected product exergy / all incoming boundary exergy: {_fmt(pct, 2)} %")
    if "fuel_performance" in report:
        fuel = report["fuel_performance"]
        pct = None if fuel["efficiency_fraction"] is None else 100.0 * fuel["efficiency_fraction"]
        lines.extend(["", "Fuel products on a lower heating value (LHV) basis",
                      f"Product LHV: {_fmt(fuel['product_lhv_kw'])} kW",
                      f"Input on declared feed and utility basis: {_fmt(fuel['energy_input_kw'])} kW",
                      f"Fuel-product energy efficiency: {_fmt(pct, 2)} %", *fuel["diagnostics"]])
    quality = report["quality"] if "quality" in report else build_quality(report)
    if quality:
        lines.extend(["", "Results to review (details and follow-up in quality.csv)",
                      f"{'Finding':<35} {'Rows':>4}   Largest absolute value or affected IDs"])
        findings = defaultdict(list)
        for row in quality:
            findings[row["issue"]].append(row)
        for issue, records in findings.items():
            values = [r for r in records if r["value"] is not None]
            if values:
                largest = max(values, key=lambda r: abs(r["value"]))
                detail = f"{largest['scope']} {largest['id']}  {_fmt(largest['value'])} {largest['unit']}"
            else:
                detail = ", ".join(row["id"] for row in records[:3])
                if len(records) > 3:
                    detail += ", ..."
            lines.append(f"{issue:<35} {len(records):>4}   {detail}")
    if plant["diagnostics"]:
        lines.extend(["", "Plant diagnostics", *plant["diagnostics"]])
    return "\n".join(lines) + "\n"


def _table(rows, preferred):
    columns = list(preferred)
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=columns)
    writer.writeheader()
    for row in rows:
        writer.writerow({k: json.dumps(v, ensure_ascii=False, allow_nan=False)
                         if isinstance(v, (dict, list, tuple)) else v for k, v in row.items()})
    return out.getvalue()


def write_results(report, destination):
    """Atomically write a ZIP with summary, JSON, and unit-labelled CSV tables."""
    path = Path(destination)
    if path.suffix.lower() != ".zip":
        raise ValueError("Results destination must end in .zip")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = None
    try:
        with NamedTemporaryFile(dir=path.parent, suffix=".tmp", delete=False) as handle:
            tmp = Path(handle.name)
            with ZipFile(handle, "w", ZIP_DEFLATED) as archive:
                if "quality" not in report:
                    # Older results can be re-exported without running a property calculation.
                    report = {**report, "quality": build_quality(report)}
                archive.writestr("summary.txt", summary(report))
                archive.writestr("analysis.json", json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False))
                if "export_summary" in report:
                    archive.writestr("export_validation.json", json.dumps(report["export_summary"], indent=2, allow_nan=False))
                archive.writestr("quality.csv", _table(report["quality"], QUALITY_COLUMNS))
                rows = []
                for s in report["streams"].values():
                    row = {k: v for k, v in s.items() if not isinstance(v, dict) and k != "substreams"}
                    if s.get("reference_option_hints"):
                        row["reference_option_hints"] = s["reference_option_hints"]
                    if "native_exergy" in s:
                        native = s["native_exergy"]
                        row.update(native_export_source=native["export_source"], native_scope=native["scope"],
                                   native_components_kw=native["components_kw"],
                                   saved_graphite_physical_kw=native["saved_graphite_physical_kw"])
                    if "export_evidence" in s:
                        evidence = s["export_evidence"]
                        row.update(aspen_whole_kw=evidence["whole_kw"], aspen_mixed_kw=evidence["mixed_kw"],
                                   aspen_solid_kw=evidence["solid_kw"], workbook_location=evidence["source_location"])
                    rows.append(row)
                archive.writestr("streams.csv", _table(rows, (
                    "id", "from_block", "to_block", "temperature_k", "pressure_pa", "mass_kg_h",
                    "enthalpy_kw", "entropy_kw_k", "physical_exergy_kw", "chemical_exergy_kw",
                    "total_exergy_kw", "split_quality", "split_method")))
                archive.writestr("blocks.csv", _table(report["blocks"], (
                    "id", "type", "energy_residual_kw", "exergy_destruction_kw",
                    "entropy_generation_kw_k", "energy_balance_status", "exergy_balance_status")))
                if report.get("groups"):
                    archive.writestr("groups.csv", _table(report["groups"].values(), (
                        "id", "members", "energy_residual_kw", "exergy_destruction_kw", "entropy_generation_kw_k")))
                archive.writestr("plant.csv", _table([report["plant"]], ("id", "energy_residual_kw",
                                                                        "exergy_destruction_kw")))
                if "fuel_performance" in report:
                    archive.writestr("fuel_performance.csv", _table([report["fuel_performance"]], (
                        "efficiency_fraction", "product_lhv_kw", "feed_lhv_kw", "energy_input_kw")))
                if "heat_recovery" in report:
                    archive.writestr("heat_recovery.csv", _table(report["heat_recovery"]["water_circuits"],
                        ("block", "role", "inlet", "outlet", "heat_gain_kw", "exergy_gain_kw")))
                if "dryer_design" in report:
                    archive.writestr("dryer_design.csv", _table([report["dryer_design"]],
                        ("outlet_ppmv", "outlet_limit_ppmv", "required_mass_per_bed_kg")))
                archive.writestr("transfers.csv", _table(report["transfers"], (
                    "block", "id", "kind", "rate_kw", "temperature_k", "exergy_kw", "basis")))
        os.replace(tmp, path)
    finally:
        if tmp is not None:
            tmp.unlink(missing_ok=True)
    return path
