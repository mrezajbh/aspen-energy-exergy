"""Public interface for Aspen archives and normalized exports."""

import json
import math
from collections import Counter
from copy import deepcopy
from importlib.metadata import PackageNotFoundError, version
from importlib.resources import files
from pathlib import Path


def _json_object(value):
    if isinstance(value, dict):
        return value
    with Path(value).open(encoding="utf-8-sig") as source:
        data = json.load(source)
    if not isinstance(data, dict):
        raise ValueError(f"Expected a JSON object in {value}")
    return data


def load_case(source):
    """Read a solved .apwz/.bkp, a normalized JSON export, or its Python dict."""
    if isinstance(source, dict) or Path(source).suffix.lower() == ".json":
        model = _json_object(source)
        if not isinstance(model.get("metadata"), dict):
            raise ValueError("Normalized data require a metadata object")
        if model["metadata"].get("enthalpy_entropy_basis", "unverified") not in ("aspen_formation", "unverified"):
            raise ValueError("Unknown H/S basis. Specify 'aspen_formation' or 'unverified'.")
    else:
        from .aspen import read_aspen, assert_cache_units
        model = read_aspen(Path(source))
        assert_cache_units(model)
    for name in ("streams", "blocks", "components"):
        rows = model.get(name)
        if not isinstance(rows, list) or not rows:
            raise ValueError(f"Case needs a nonempty {name} list")
        if any(not isinstance(r, dict) or not isinstance(r.get("id"), str) or not r["id"] for r in rows):
            raise ValueError(f"Each entry in {name} must have a nonempty string id")
        if len({r["id"] for r in rows}) != len(rows):
            raise ValueError(f"Duplicate id in {name}")
    return model


def _settings(value, model=None):
    if value is None or value == "as-saved":
        return {}
    if value == "biogas-ft":
        profile = files("aspen_energy_exergy").joinpath("profiles", "biogas-ft.json")
        settings = json.loads(profile.read_text(encoding="utf-8"))
        from .aspen import EXERGY_REPORT_BACKUP
        if model and model["metadata"].get("backup_sha256") == EXERGY_REPORT_BACKUP:
            _refresh_report_profile(settings, model)
        return settings
    return _json_object(value)


def _refresh_report_profile(settings, model):
    """Apply the audited utility/power conventions to the revised operating states."""
    settings["backup_sha256"] = model["metadata"]["backup_sha256"]
    settings["name"] = "Biogas to Jetfuel - FT | November exergy report and declared utility boundaries"
    streams = {s["id"]: s for s in model["streams"]}
    blocks = {b["id"]: b for b in model["blocks"]}
    donor = blocks["H-REC-HX"]
    donor_heat = next(t["rate_kw"] for t in donor["transfers"] if t["id"] == "QQ")
    donor_temperature = donor_heat / (streams["COLDER2"]["entropy_kw_k"] - streams["HCK-OUT"]["entropy_kw_k"])
    for contact in settings["transfers"]:
        bid, tid = contact["block"], contact["id"]
        if contact["basis"] in {"Uniform process outlet temperature proxy", "Flue gas outlet temperature proxy"}:
            values = [streams[s]["temperature_k"] for s in blocks[bid]["outlets"] if streams[s]["mass_kg_h"]]
            if max(values) - min(values) > 1e-6:
                raise ValueError(f"{bid}: multiple outlet temperatures for a uniform heat boundary")
            contact["temperature_k"] = values[0]
        elif bid == "COOLER":
            contact["temperature_range_k"] = [streams[s]["temperature_k"] for s in ("WAT", "WAT-C")]
        elif tid == "COOL-ST":
            contact["temperature_k"] = streams["LPS-HEAD"]["temperature_k"]
        elif tid == "QQ":
            contact["temperature_k"] = donor_temperature
        elif bid == "T-601":
            contact["temperature_k"] = next(t["process_temperature_k"] for t in blocks[bid]["transfers"]
                                             if t["id"] == tid)


def _extend_profile(settings, model=None):
    """Resolve a named profile while retaining each unmodified transfer field."""
    name = settings["base_profile"]
    if name != "biogas-ft":
        raise ValueError(f"Unknown base_profile: {name!r}")
    base = _settings(name, model)
    if settings.get("backup_sha256", base["backup_sha256"]) != base["backup_sha256"]:
        raise ValueError("Cannot change the backup SHA-256 of an inherited case profile")
    overrides = settings.get("transfers", [])
    if not isinstance(overrides, list):
        raise ValueError("settings.transfers must be a list")
    index = {(row["block"], row["id"]): row for row in base["transfers"]}
    seen = set()
    for row in overrides:
        if not isinstance(row, dict) or not all(isinstance(row.get(k), str) for k in ("block", "id")):
            raise ValueError("Each transfer override needs a block and transfer id")
        key = (row["block"], row["id"])
        if key in seen:
            raise ValueError(f"Duplicate transfer override: {key}")
        seen.add(key)
        inherited = dict(index.get(key, {}))
        if not isinstance(row.get("basis"), str) or not row["basis"].strip():
            raise ValueError(f"Identify the source or assumption for {key} in 'basis'")
        if row.get("kind") == "work" or "temperature_range_k" in row:
            inherited.pop("temperature_k", None)
        if row.get("kind") == "work" or "temperature_k" in row:
            inherited.pop("temperature_range_k", None)
        index[key] = {**inherited, **row}
    merged = {**base, **settings, "transfers": list(index.values())}
    del merged["base_profile"]
    return merged


def _apply_settings(model, settings):
    unknown = settings.keys() - {"name", "backup_sha256", "transfers", "product_streams",
                                 "confirm_formation_datum", "fuel_feeds", "utility_pairs", "groups"}
    if unknown:
        raise ValueError(f"Unknown settings: {sorted(unknown)}")
    if "backup_sha256" in settings and settings["backup_sha256"] != model["metadata"].get("backup_sha256"):
        raise ValueError("The settings belong to a different Aspen backup (SHA-256 mismatch)")
    confirmation = settings.get("confirm_formation_datum", False)
    if type(confirmation) is not bool:
        raise ValueError("confirm_formation_datum must be a boolean")
    if confirmation:
        model = {**model, "metadata": {**model["metadata"], "enthalpy_entropy_basis": "aspen_formation",
                                       "datum_confirmation": "User-confirmed in analysis settings"}}
    overrides = settings.get("transfers", [])
    if not isinstance(overrides, list):
        raise ValueError("settings.transfers must be a list")
    if not overrides:
        return model
    blocks = deepcopy(model["blocks"])
    targets = {(b["id"], t["id"]): t for b in blocks for t in b["transfers"]}
    seen = set()
    for entry in overrides:
        if not isinstance(entry, dict) or entry.keys() - {"block", "id", "kind", "temperature_k",
                                                         "temperature_range_k", "basis"}:
            raise ValueError(f"Invalid transfer override: {entry}")
        key = (entry.get("block"), entry.get("id"))
        if key not in targets or key in seen:
            raise ValueError(f"Unknown or duplicate transfer override: {key}")
        seen.add(key)
        if not isinstance(entry.get("basis"), str) or not entry["basis"].strip():
            raise ValueError(f"Identify the source or assumption for {key} in 'basis'")
        transfer = targets[key]
        if "kind" in entry:
            if entry["kind"] not in ("work", "heat"):
                raise ValueError(f"Invalid transfer kind at {key}")
            transfer["kind"] = entry["kind"]
        if transfer["kind"] == "work":
            transfer.pop("temperature_k", None)
            transfer.pop("temperature_range_k", None)
        if "temperature_range_k" in entry:
            if "temperature_k" in entry or transfer["kind"] != "heat" or entry["temperature_range_k"] is None:
                raise ValueError(f"Specify a temperature range only on heat with no scalar temperature at {key}")
            transfer.pop("temperature_k", None)
            transfer["temperature_range_k"] = entry["temperature_range_k"]
        if "temperature_k" in entry:
            temp = entry["temperature_k"]
            if (transfer["kind"] != "heat" or type(temp) not in (float, int)
                    or not math.isfinite(temp) or temp <= 0):
                raise ValueError(f"Expected a positive heat boundary temperature in K at {key}")
            transfer["temperature_k"] = temp
            transfer.pop("temperature_range_k", None)
        transfer["basis"] = entry["basis"]
    return {**model, "blocks": blocks}


def write_native_exergy_template(source, destination):
    """Create a case-bound CSV for returning Aspen EXERGYFL values in kW."""
    from .native_exergy import write_native_exergy_template as write_template

    return write_template(load_case(source), destination)


def _analysis_inputs(source, settings, references, native_exergy, report_text=None):
    model = load_case(source)
    rates = {}
    workbook = native_exergy is not None and Path(native_exergy).suffix.lower() == ".xlsx"
    if report_text is not None and not workbook:
        raise ValueError("Use the text report with a native .xlsx workbook")
    if workbook:
        from .aspen_exergy_report import parse_export
        if isinstance(source, dict) or Path(source).suffix.lower() not in (".apwz", ".bkp"):
            raise ValueError("A native Aspen workbook needs its associated .apwz or .bkp to verify operating states")
        export = parse_export(native_exergy, report_text, source, case_model=model)
        by_id = export["streams"]
        if set(by_id) != {s["id"] for s in model["streams"]}:
            raise ValueError("Workbook and archive stream identifiers differ")
        rows = []
        for saved in model["streams"]:
            evidence = by_id[saved["id"]]
            values = evidence["operating"]
            # The detailed workbook includes empty substream sections on every column.
            # Keep only active substreams in the normalized case, as in the archive reader.
            active_subs = {k: sub for k, sub in values["substreams"].items() if sub["molar_kmol_h"] > 0}
            row = {**saved, **values, "substreams": active_subs,
                   "export_evidence": {k: v for k, v in evidence.items() if k != "operating"}}
            rows.append(row)
            if values["mass_kg_h"] == 0:
                if evidence["whole_kw"] not in (None, 0.0):
                    raise ValueError(f"{saved['id']}: a zero-flow stream has nonzero EXERGYFL")
                continue
            parts = {label: evidence[key] for label, key in (("MIXED", "mixed_kw"), ("CISOLID", "solid_kw"))
                     if evidence[key] is not None}
            rates[saved["id"]] = dict(physical_exergy_kw=evidence["whole_kw"], components_kw=parts,
                scope="whole_stream", saved_graphite_physical_kw=None,
                export_source=f"{Path(native_exergy).name}, {evidence['source_location']}",
                case_sha256=model["metadata"]["backup_sha256"], reference_temperature_k=298.15,
                reference_pressure_pa=101325.0)
        metadata = {**model["metadata"], "source_properties": "Aspen workbook material states and native EXERGYFL",
                    "material_state_source": Path(native_exergy).name,
                    "equipment_source": model["metadata"]["source_filename"]}
        model = {**model, "streams": rows, "export_summary": export["summary"], "metadata": metadata}
    selection = "explicit" if settings is not None else "as saved"
    if settings is None and not isinstance(source, dict) and Path(source).suffix.lower() in {".apwz", ".bkp"}:
        known = _settings("biogas-ft", model)
        if model["metadata"].get("backup_sha256") == known["backup_sha256"]:
            settings = known
            selection = "biogas-ft matched by native backup SHA-256"
    settings = _settings(settings, model)
    if "base_profile" in settings:
        selection = f"explicit overrides of {settings['base_profile']}"
        settings = _extend_profile(settings, model)
    model = _apply_settings(model, settings)
    refs = []
    if references is None:
        references = []
    elif not isinstance(references, (list, tuple)):
        references = [references]
    for item in references:
        if not isinstance(item, dict) and Path(item).suffix.lower() == ".csv":
            from .flash_exchange import read_flash_returns
            refs.append(read_flash_returns(model, item))
        else:
            refs.append(_json_object(item))
    if native_exergy is not None and not workbook:
        from .native_exergy import read_native_exergy
        rates = read_native_exergy(native_exergy, model)
    return model, settings, selection, refs, rates


def plan_reference_flashes(source, *, settings=None, references=None, native_exergy=None, report_text=None):
    """Plan outstanding native reference flashes from saved states without running a property model."""
    from .reference_plan import build_reference_plan

    model, _, _, refs, rates = _analysis_inputs(source, settings, references, native_exergy, report_text)
    return build_reference_plan(model, references=refs, native_exergy=rates)


def analyze(source, *, settings=None, references=None, estimate_splits=True, native_exergy=None, report_text=None,
            balance_engine="builtin"):
    """Return stream, equipment, subsystem and plant balances on the 298.15 K, 101325 Pa basis.

    references accepts matching ambient-flash dictionaries, JSON documents or flash-return CSVs.
    Pass a list for nonoverlapping partial returns from several batches.
    native_exergy accepts an Aspen EXERGYFL CSV or a full material report workbook;
    report_text adds the printed Aspen report comparison for the workbook. The supplied backup selects
    the biogas-ft settings automatically. Use settings='as-saved' for the raw transfer inventory.
    balance_engine='exerpy-tespy' evaluates and publishes the upstream imported-state
    balances; it requires estimate_splits=False and never runs a property flash.
    The default builtin engine retains the historical API behavior.
    """
    from .balances import calculate_balances
    from .properties import ELEMENT_EXERGY_J_MOL, P0_PA, T0_K, calculate_streams
    from .quality import build_quality

    if type(estimate_splits) is not bool:
        raise ValueError("estimate_splits must be a boolean")
    if not isinstance(balance_engine, str) or balance_engine not in {"builtin", "exerpy-tespy"}:
        raise ValueError("balance_engine must be builtin or exerpy-tespy")
    if balance_engine == "exerpy-tespy" and estimate_splits:
        raise ValueError("ExerPy/TESPy imported-state mode requires estimate_splits=False; supply Aspen states.")
    model, settings, selection, refs, rates = _analysis_inputs(source, settings, references, native_exergy, report_text)
    properties = calculate_streams(model, references=refs, estimate_splits=estimate_splits, native_exergy=rates)
    result = calculate_balances(model, properties, product_streams=settings.get("product_streams"),
                                groups=settings.get("groups"))
    if "fuel_feeds" in settings:
        from .performance import fuel_performance
        result["fuel_performance"] = fuel_performance(model, result["plant"], settings["fuel_feeds"],
                                                      settings.get("utility_pairs", ()))
    elif "utility_pairs" in settings:
        raise ValueError("Specify fuel_feeds when selecting material utility_pairs for fuel performance")
    streams = {s["id"]: {**s, **properties[s["id"]]} for s in model["streams"]}
    versions = {}
    for package in ("aspen-energy-exergy", "thermo", "chemicals", "fluids", "iapws"):
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            continue
    report = {
        "schema_version": 1,
        "source": model["metadata"],
        "environment": {"temperature_k": T0_K, "pressure_pa": P0_PA,
                        "element_exergy_j_per_mol_atoms": dict(ELEMENT_EXERGY_J_MOL),
                        "reference": "Tang et al. (2021), doi:10.3303/CET2183089, equation 5 and table 1",
                        "enthalpy_entropy_basis": model["metadata"].get("enthalpy_entropy_basis", "unverified")},
        "software_versions": versions,
        "settings": settings,
        "settings_selection": selection,
        "estimate_splits": estimate_splits,
        "native_exergy_import": None if native_exergy is None else {
            "file": Path(native_exergy).name, "accepted_streams": len(rates)
        },
        "components": model["components"],
        "streams": streams,
        "coverage": {
            "streams": len(streams),
            "total_exergy_available": sum(s["total_exergy_kw"] is not None for s in streams.values()),
            "splits_available": sum(s["physical_exergy_kw"] is not None and s["chemical_exergy_kw"] is not None
                                    for s in streams.values()),
            "split_quality": dict(Counter(s["split_quality"] for s in streams.values())),
            "blocks": len(result["blocks"]),
            "block_exergy_balances_available": sum(b["exergy_destruction_kw"] is not None for b in result["blocks"]),
        },
        "source_warnings": model.get("warnings", []),
        "energy_links": model.get("energy_links", []),
        **result,
    }
    if "export_summary" in model:
        report["export_summary"] = model["export_summary"]
        report["environment"]["physical_reference_basis"] = model["export_summary"]["reference_basis"]
    if balance_engine == "exerpy-tespy":
        from .engines import promote_external_results, run_external_engines
        report["external_engines"] = run_external_engines(model, properties, result)
        promote_external_results(report)
    report["balance_engine"] = balance_engine
    report["quality"] = build_quality(report)
    return report
