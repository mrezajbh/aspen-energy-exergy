"""Exchange physical exergy flow rates with an Aspen EXERGYFL stream report."""

import csv
import io
from math import fsum, isclose, isfinite
from pathlib import Path

from ._estimate import is_graphite
from ._exchange import case_hashes as _case_hashes, digest as _digest, publish_template
from .properties import P0_PA, T0_K


_RATES = {"whole_stream": "whole_stream_EXERGYFL_kW",
          "MIXED": "MIXED_EXERGYFL_kW", "CISOLID": "CISOLID_EXERGYFL_kW"}
_STATE = ("operating_T_C", "operating_P_bar", "operating_mass_kg_h", "operating_H_kW", "operating_S_kW_K")
_HEADERS = ("id", "property_method", "from_block", "to_block", *_STATE, "required_substreams",
            "saved_graphite_physical_kW", "target_exergy_reference_T_K", "target_exergy_reference_P_Pa",
            "original_case_sha256", "export_source", "same_operating_state_confirmed", "reference_T_K",
            "reference_P_Pa", "scope", *_RATES.values())


def _saved_graphite(stream, components, formation):
    solid = (stream.get("substreams") or {}).get("CISOLID", {})
    active = [cid for cid, n in solid.get("component_flows_kmol_h", {}).items() if n > 0]
    if not formation or len(active) != 1 or not is_graphite(components[active[0]]):
        return None
    h, s = solid.get("enthalpy_kw"), solid.get("entropy_kw_k")
    return None if h is None or s is None else h - T0_K * s


def _source_fields(stream, graphite):
    t, p = stream.get("temperature_k"), stream.get("pressure_pa")
    subs = stream.get("substreams") or {}
    # Substream rates are usable only when the export identifies every active
    # substream. Whole-stream rates also work for unknown or new substreams.
    known_subs = bool(subs) and set(subs) <= {"MIXED", "CISOLID"}
    required = " + ".join(k for k in ("MIXED", "CISOLID") if k in subs) if known_subs else "whole_stream"
    return {
        "id": stream["id"], "property_method": stream.get("property_method") or "",
        "from_block": stream.get("from_block") or "", "to_block": stream.get("to_block") or "",
        "operating_T_C": t - 273.15 if t is not None else None,
        "operating_P_bar": p / 1e5 if p is not None else None,
        "operating_mass_kg_h": stream.get("mass_kg_h"), "operating_H_kW": stream.get("enthalpy_kw"),
        "operating_S_kW_K": stream.get("entropy_kw_k"), "required_substreams": required,
        "saved_graphite_physical_kW": graphite, "target_exergy_reference_T_K": T0_K,
        "target_exergy_reference_P_Pa": P0_PA,
    }


def _number(value, field):
    try:
        num = float(value)
    except (ValueError, TypeError):
        raise ValueError(f"{field} requires a numeric value") from None
    if not isfinite(num):
        raise ValueError(f"{field} requires a finite value")
    return num


def write_native_exergy_template(model, destination):
    """Write a blank, case-bound return sheet for streams needing a physical split.

    model is a normalized Aspen case. Source state and units are prefilled.
    An existing identical template is reused; a filled sheet is preserved.
    """
    path = Path(destination)
    if path.suffix.lower() != ".csv":
        raise ValueError("Use a .csv destination for the EXERGYFL template")
    digest = _case_hashes(model)[0]
    components = {c["id"]: c for c in model["components"]}
    formation = model["metadata"].get("enthalpy_entropy_basis") == "aspen_formation"
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=(*_HEADERS, "operating_state_sha256"))
    writer.writeheader()
    for stream in model["streams"]:
        if stream.get("molar_kmol_h") == 0:
            continue
        active = [cid for cid, n in (stream.get("component_flows_kmol_h") or {}).items() if n > 0]
        reference_graphite = (formation and len(active) == 1 and is_graphite(components[active[0]])
                              and all(stream.get(k) is not None for k in ("enthalpy_kw", "entropy_kw_k")))
        if reference_graphite:
            continue
        row = _source_fields(stream, _saved_graphite(stream, components, formation))
        writer.writerow({**row, "original_case_sha256": digest, "operating_state_sha256": _digest(stream)})
    data = output.getvalue().encode("utf-8-sig")
    return publish_template(path, data)


def _check_source(row, expected, stream, hashes, native_backup_hash):
    sid = stream["id"]
    if row["original_case_sha256"].lower() not in hashes:
        raise ValueError(f"{sid}: EXERGYFL sheet belongs to a different case")
    if "operating_state_sha256" in row:
        fingerprint = row["operating_state_sha256"].lower()
        if fingerprint != _digest(stream):
            # v0.2 native templates predate extraction of producer option hints.
            # The backup hash still binds the actual producer options in the deck.
            legacy = {k: v for k, v in stream.items() if k not in ("reference_option_hints", "reference_options_key")}
            if row["original_case_sha256"].lower() != native_backup_hash or fingerprint != _digest(legacy):
                raise ValueError(f"{sid}: operating state fingerprint differs")
    for key in ("property_method", "from_block", "to_block", "required_substreams"):
        if row[key] != expected[key]:
            raise ValueError(f"{sid}: {key} differs from the saved operating state")
    for field in (*_STATE, "target_exergy_reference_T_K", "target_exergy_reference_P_Pa",
                  "saved_graphite_physical_kW"):
        value = expected[field]
        if field == "saved_graphite_physical_kW" and not row[field]:
            # A template prepared before the formation datum was confirmed
            # can omit this informational hint. Always compute it from source.
            continue
        if value is None and not row[field]:
            continue
        if value is None or not isclose(_number(row[field], f"{sid} {field}"), value, rel_tol=1e-10, abs_tol=1e-8):
            raise ValueError(f"{sid}: {field} differs from the case")


def _physical_rate(row, stream, graphite):
    sid, scope = stream["id"], row["scope"]
    rates = {kind: _number(row[key], f"{sid} {key}") for kind, key in _RATES.items() if row[key]}
    active = set(stream.get("substreams") or {})
    used_graphite = None
    if scope == "whole_stream":
        if rates.keys() != {"whole_stream"}:
            raise ValueError(f"{sid}: whole_stream needs exactly one whole-material EXERGYFL rate in kW")
    elif scope == "mixed_with_saved_graphite":
        if rates.keys() != {"MIXED"} or active != {"MIXED", "CISOLID"} or graphite is None:
            raise ValueError(f"{sid}: mixed_with_saved_graphite needs a MIXED rate and verified native graphite")
        used_graphite = graphite
    elif scope == "substreams":
        if not active or not active <= {"MIXED", "CISOLID"}:
            raise ValueError(f"{sid}: substream inventory is unavailable or unsupported; use whole_stream")
        # Earlier MIXED-only returns permit an explicitly entered zero for the
        # absent CISOLID. A nonzero value would include an unrelated stream.
        if "CISOLID" not in active and rates.get("CISOLID") == 0:
            del rates["CISOLID"]
        if rates.keys() != active:
            raise ValueError(f"{sid}: enter one EXERGYFL kW rate for each active substream {sorted(active)}")
    else:
        raise ValueError(f"{sid}: select whole_stream, substreams, or mixed_with_saved_graphite")
    try:
        rate = fsum([*rates.values(), used_graphite or 0.0])
    except OverflowError:
        raise ValueError(f"{sid}: summed physical exergy is outside the finite numeric range") from None
    if stream.get("molar_kmol_h") == 0 and rate != 0:
        raise ValueError(f"{sid}: nonzero EXERGYFL on a zero-flow stream")
    return rate, rates, used_graphite


def read_native_exergy(path, model):
    """Read completed EXERGYFL kW rows. Unfilled rows keep their existing methods.

    Both the original project return sheet and library-generated templates are
    accepted. Checks apply to every row, including rows without an Aspen value.
    """
    streams = {s["id"]: s for s in model["streams"]}
    components = {c["id"]: c for c in model["components"]}
    formation = model["metadata"].get("enthalpy_entropy_basis") == "aspen_formation"
    hashes, seen, result = _case_hashes(model), set(), {}
    with Path(path).open(encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file, strict=True)
        try:
            names = reader.fieldnames or []
            if len(names) != len(set(names)) or not set(_HEADERS) <= set(names):
                raise ValueError("Invalid EXERGYFL CSV headers. Use the return template with EXERGYFL rates in kW.")
            for values in reader:
                if None in values or None in values.values():
                    raise ValueError(f"Malformed EXERGYFL row at line {reader.line_num}")
                row = {k: v.strip() for k, v in values.items()}
                if not any(row.values()):
                    continue
                sid = row["id"]
                if sid in seen or sid not in streams:
                    raise ValueError(f"Duplicate or unknown EXERGYFL stream ID {sid!r} at line {reader.line_num}")
                seen.add(sid)
                stream = streams[sid]
                graphite = _saved_graphite(stream, components, formation)
                _check_source(row, _source_fields(stream, graphite), stream, hashes,
                              model["metadata"].get("backup_sha256") if model["metadata"].get("source_properties")
                              == "Aspen saved ADS SI bulk properties" else None)
                if not row["scope"] and not any(row[k] for k in _RATES.values()):
                    continue
                if not row["export_source"] or row["same_operating_state_confirmed"].lower() != "yes":
                    raise ValueError(f"{sid}: provide export_source and confirm the same operating state with yes")
                for key, target, tol in (("reference_T_K", T0_K, 1e-6), ("reference_P_Pa", P0_PA, 1e-4)):
                    if not isclose(_number(row[key], f"{sid} {key}"), target, abs_tol=tol, rel_tol=0):
                        raise ValueError(f"{sid}: {key} differs from the analysis reference environment")
                rate, parts, used_graphite = _physical_rate(row, stream, graphite)
                result[sid] = dict(physical_exergy_kw=rate, components_kw=parts, saved_graphite_physical_kw=used_graphite,
                                   scope=row["scope"], export_source=row["export_source"],
                                   case_sha256=row["original_case_sha256"].lower(),
                                   reference_temperature_k=T0_K, reference_pressure_pa=P0_PA)
        except csv.Error as error:
            raise ValueError(f"Invalid EXERGYFL CSV near line {reader.line_num}: {error}") from error
    return result
