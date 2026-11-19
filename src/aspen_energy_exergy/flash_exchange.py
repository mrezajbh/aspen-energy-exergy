"""Portable Aspen flash requests and compact molar-property returns."""

import csv
import io
import json
from math import isfinite
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

from ._exchange import case_hashes, digest, publish_template
from .properties import REFERENCE_PLAN_SCHEMA
from .report import _table


_IDENTITY = ("id", "applies_to", "property_method", "case_sha256", "operating_state_sha256")
_MEASUREMENTS = ("aspen_reference_id", "source_location", "T_K", "P_Pa", "h0_J_mol", "s0_J_molK",
                 "representative_verified", "peers_verified")

_INSTRUCTIONS = """Aspen Plus ambient flash requests

This packet requests equilibrium reference states at 298.15 K and 101325 Pa.
Each row in flash_returns.csv is a candidate flash. flash_groups.csv connects
it to the operating streams. flash_compositions.csv lists the nonzero component
rates for each representative in kmol/h, separated by substream when available.

1. Flash a copy of the representative composition at the reference conditions.
   Retain the Aspen property method, component and interaction data, substream
   allocations and phase options, including free-water settings. The packet
   identifies saved option declarations. Check all settings in the simulation.
2. In flash_returns.csv fill in aspen_reference_id, source_location, the actual
   T_K and P_Pa, and whole-material molar enthalpy and entropy in h0_J_mol and
   s0_J_molK. Aspen kJ/kmol and kJ/(kmol K) have the same numerical values
   as J/mol and J/(mol K). For mixed fluid/solid streams include both parts.
   Use the same enthalpy and entropy convention as the operating case.
3. Enter yes in representative_verified after confirming composition, substream
   allocation, property settings, equilibrium phases and units in the flash.
   In applies_to retain only peers whose reference settings and composition you
   have verified. For any retained peers enter yes in peers_verified. To use
   a flash only for the representative enter [] in applies_to and leave
   peers_verified blank. Return fields for outstanding flashes can stay blank.
4. Supply the filled CSV to the analysis with --references flash_returns.csv.
   Use --native-splits-only to see the rows supported by reference data.

flash_plan.json is an alternative return format. Its default flow_SI fields
are kW and kW/K at the listed component rates. Fill one return format.
"""


def write_flash_plan(plan, destination):
    """Write a fillable JSON plan or a ZIP containing the human handoff and CSVs."""
    if not isinstance(plan, dict) or plan.get("plan_schema") != REFERENCE_PLAN_SCHEMA:
        raise ValueError("Expected a reference flash plan from plan_reference_flashes()")
    path = Path(destination)
    document = (json.dumps(plan, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    if path.suffix.lower() == ".json":
        return publish_template(path, document)
    if path.suffix.lower() != ".zip":
        raise ValueError("Use .zip for a flash request packet or .json for the fillable plan")
    returns, groups, components = [], [], []
    for record in plan["streams"]:
        sid = record["id"]
        returns.append({"id": sid, "applies_to": record["applies_to"], "property_method": record["property_method"],
                        "case_sha256": plan["case_sha256"], "operating_state_sha256": record["operating_state_sha256"]})
        groups.append({k: record[k] for k in ("id", "member_count", "member_ids", "property_method",
                                              "operating_phase_sets", "reference_option_hints", "member_sources", "notes")})
        substreams = record["substream_component_flows_kmol_h"]
        if not substreams or any(rates is None for rates in substreams.values()):
            substreams = {"whole_stream": dict(zip(record["component_ids"], record["component_flows_kmol_h"]))}
        for sub, rates in substreams.items():
            for component, rate in rates.items():
                if rate > 0:
                    components.append({"id": sid, "substream": sub, "component": component, "flow_kmol_h": rate})
    contents = {
        "instructions.txt": (_INSTRUCTIONS + f"\nSource case: {plan.get('source_filename') or 'normalized Aspen export'}"
                             f"\nCase SHA-256: {plan['case_sha256']}\n"),
        "flash_plan.json": document,
        "flash_groups.csv": _table(groups, ("id", "member_count", "member_ids", "property_method")),
        "flash_compositions.csv": _table(components, ("id", "substream", "component", "flow_kmol_h")),
        "flash_returns.csv": _table(returns, (*_IDENTITY, *_MEASUREMENTS)),
    }
    buffer = io.BytesIO()
    with ZipFile(buffer, "w") as archive:
        for filename, content in contents.items():
            # Fixed ZIP timestamps make repeated identical requests idempotent.
            archive.writestr(ZipInfo(filename), content, compress_type=ZIP_DEFLATED)
    return publish_template(path, buffer.getvalue())


def read_flash_returns(model, path):
    """Bind a partial CSV to the source and convert it to a validated reference document."""
    hashes = case_hashes(model)
    streams = {s["id"]: s for s in model["streams"]}
    component_ids = [c["id"] for c in model["components"]]
    records = []
    seen = set()
    with Path(path).open(encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file, strict=True)
        try:
            fields = reader.fieldnames or []
            if len(set(fields)) != len(fields) or not set((*_IDENTITY, *_MEASUREMENTS)) <= set(fields):
                raise ValueError("Expected the flash_returns.csv header from a reference flash packet")
            for row in reader:
                if None in row or None in row.values():
                    raise ValueError(f"Malformed flash return at CSV line {reader.line_num}")
                row = {k: v.strip() for k, v in row.items()}
                if not any(row.values()):
                    continue
                sid = row["id"]
                if sid not in streams or sid in seen:
                    raise ValueError(f"Unknown or duplicate representative {sid!r}")
                seen.add(sid)
                stream = streams[sid]
                if row["case_sha256"].lower() not in hashes or row["operating_state_sha256"].lower() != digest(stream):
                    raise ValueError(f"{sid}: flash return does not match this case and operating state")
                if row["property_method"] != stream.get("property_method"):
                    raise ValueError(f"{sid}: flash return property method differs")
                if not any(row[k] for k in _MEASUREMENTS):
                    continue
                peers = json.loads(row["applies_to"] or "[]")
                if not isinstance(peers, list) or not all(isinstance(peer, str) for peer in peers):
                    raise ValueError(f"{sid}: applies_to must be a JSON array of stream IDs")
                if row["representative_verified"].lower() != "yes":
                    raise ValueError(f"{sid}: confirm the flash conditions and composition in representative_verified")
                if peers and row["peers_verified"].lower() != "yes":
                    raise ValueError(f"{sid}: confirm shared settings in peers_verified or use [] for applies_to")
                values = {}
                for key in ("T_K", "P_Pa", "h0_J_mol", "s0_J_molK"):
                    try:
                        values[key] = float(row[key])
                    except ValueError as exc:
                        raise ValueError(f"{sid}: {key} must contain the returned number") from exc
                    if not isfinite(values[key]):
                        raise ValueError(f"{sid}: {key} must be finite")
                records.append({
                    "id": sid, "applies_to": peers, "property_method": row["property_method"],
                    "aspen_reference_id": row["aspen_reference_id"], "source_location": row["source_location"],
                    "component_ids": component_ids, "component_flows_kmol_h": [
                        stream["component_flows_kmol_h"].get(k, 0.) for k in component_ids],
                    "quantity_basis": "molar_SI", "phase_convention": "equilibrium",
                    "same_aspen_property_options": True, "same_aspen_property_options_for_all": bool(peers),
                    **values,
                })
        except csv.Error as exc:
            raise ValueError(f"Invalid flash return CSV near line {reader.line_num}: {exc}") from exc
    return {"case_sha256": hashes[0], "streams": records}
