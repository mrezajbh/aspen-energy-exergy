"""Trace calculation and source-data findings back to an equipment item or stream."""

from math import isclose

from .properties import _reference_profile, _same_fractions

QUALITY_COLUMNS = ("issue", "scope", "id", "metric", "value", "unit", "evidence", "action")


def _finding(issue, scope, item, metric, value, unit, evidence, action):
    return dict(issue=issue, scope=scope, id=item, metric=metric, value=value, unit=unit,
                evidence=evidence, action=action)


def _boundary_findings(scope, boundary, groups):
    bid = boundary["id"]
    diagnostics = boundary.get("diagnostics", [])
    if boundary.get("exergy_balance_status") == "negative":
        evidence = "The signed exergy balance is below the -1 W reporting threshold."
        if scope == "block":
            enclosing = [g for g in groups if bid in g["members"] and len(g["members"]) > 1
                         and g.get("exergy_destruction_kw") is not None]
            if enclosing:
                g = min(enclosing, key=lambda g: (len(g["members"]), g["id"]))
                evidence += f" Enclosing group {g['id']} has signed destruction {g['exergy_destruction_kw']:+.3f} kW."
        yield _finding("Negative exergy destruction", scope, bid, "exergy_destruction_kw",
                       boundary["exergy_destruction_kw"], "kW", evidence,
                       "Check the Aspen inlet and outlet states, signs on heat and work, and heat-contact temperatures.")
    elif boundary.get("t0_entropy_generation_kw") is not None and boundary["t0_entropy_generation_kw"] < -0.001:
        yield _finding("Negative entropy generation", scope, bid, "t0_entropy_generation_kw",
                       boundary["t0_entropy_generation_kw"], "kW", "T0 times entropy generation is below -1 W.",
                       "Check the saved entropy rates and temperature of each heat transfer.")
    if boundary.get("energy_balance_status") == "residual":
        tol = boundary["energy_balance_tolerance_kw"]
        yield _finding("First-law residual", scope, bid, "energy_residual_kw", boundary["energy_residual_kw"],
                       "kW", f"The absolute residual exceeds the boundary tolerance of {tol:.6g} kW.",
                       "Check the material enthalpy rates and every heat and work transfer across this boundary.")
    if "Material mass inflow and outflow differ" in diagnostics:
        yield _finding("Mass imbalance", scope, bid, "mass_residual_kg_h", boundary["mass_residual_kg_h"],
                       "kg/h", "The material flow balance exceeds the mass tolerance.",
                       "Verify crossing stream rates, ports and active substreams.")
    if boundary.get("element_balance_status") == "residual":
        residuals = boundary["element_residual_kmol_h"]
        element = max(residuals, key=lambda e: abs(residuals[e]))
        yield _finding("Elemental imbalance", scope, bid, f"{element} residual", residuals[element],
                       "kmol atoms/h", "The boundary failed elemental closure. The value is its largest absolute residual, with sign.",
                       "Check component formulas, component flow rates and reaction stoichiometry.")
    if any(message.startswith("Exergy balance differs from T0 Sgen") for message in diagnostics):
        yield _finding("Exergy identity residual", scope, bid, "exergy_identity_residual_kw",
                       boundary["exergy_identity_residual_kw"], "kW",
                       "Exergy destruction differs from the sum of lost work, energy residual and elemental-reference residual.",
                       "Verify that stream properties use a common formation datum and include all substreams.")


def _stream_findings(stream):
    sid = stream["id"]
    if stream.get("native_exergy") is not None and stream["physical_exergy_kw"] < 0:
        yield _finding("Signed negative Aspen EXERGYFL", "stream", sid, "physical_exergy_kw",
                       stream["physical_exergy_kw"], "kW", stream["native_exergy"]["export_source"],
                       "Retain Aspen's exported sign when subtracting physical exergy from total exergy.")
    if stream.get("reference_difference_kw") is not None:
        yield _finding("Conflicting exergy references", "stream", sid, "reference_difference_kw",
                       stream["reference_difference_kw"], "kW",
                       f"The chemical rate differs from reference stream {stream['comparison_reference_id']} "
                       "after scaling to the same throughput and composition.",
                       "Check the substream scope, reference environment and property options of both exports.")
    missing = [k for k in ("total_exergy_kw", "physical_exergy_kw", "chemical_exergy_kw") if stream.get(k) is None]
    if missing:
        notes = " ".join(stream.get("notes", []))
        action = ("Verify the formation basis, elemental composition, and whole-stream enthalpy and entropy rates."
                  if "total_exergy_kw" in missing else
                  "Supply Aspen EXERGYFL for this stream or a matching same-composition ambient flash.")
        yield _finding("Incomplete stream exergy", "stream", sid, ", ".join(missing), None, "kW",
                       notes or "The operating record lacks a required input.", action)
    if stream.get("split_quality") != "estimated":
        return
    method = stream.get("independent_split_method") or stream.get("split_method", "")
    if "wilson_k" in method:
        yield _finding("Wilson phase approximation", "stream", sid, "independent_split_method", None, "",
                       f"The Peng–Robinson equilibrium flash used the {method} fallback.",
                       "Compare the phase split and physical exergy with an Aspen flash at the same composition.")
    phases = set(stream.get("active_phases") or ())
    if {"Liquid1", "Liquid2"}.issubset(phases) and method.startswith("thermo_"):
        yield _finding("Two-liquid phase mismatch", "stream", sid, "saved_liquid_phases", 2, "phases",
                       "Aspen saved two liquid phases. The independent stream estimate allows one liquid phase.",
                       "Use a matching Aspen reference with both liquid phases enabled.")
    if any("heat-capacity correlation is extrapolated" in note for note in stream.get("notes", [])):
        yield _finding("Heat-capacity extrapolation", "stream", sid, "operating_temperature_k",
                       stream.get("temperature_k"), "K", "At least one heat-capacity correlation was extrapolated between the operating and ambient states.",
                       "Check property validity over the temperature range or supply Aspen EXERGYFL.")


def _cross_option_comparison(block, streams):
    inlets, outlets = block.get("inlets", ()), block.get("outlets", ())
    if len(inlets) != 1 or len(outlets) != 1:
        return None
    first, last = streams[inlets[0]], streams[outlets[0]]
    keys = first.get("reference_options_key"), last.get("reference_options_key")
    method = first.get("property_method")
    if not all(keys) or keys[0] == keys[1] or not method or last.get("property_method") != method:
        return None
    if not _same_fractions(_reference_profile(first), _reference_profile(last)):
        return None
    n_in, n_out = first["molar_kmol_h"], last["molar_kmol_h"]
    e_in, e_out = first.get("chemical_exergy_kw"), last.get("chemical_exergy_kw")
    if e_in is None or e_out is None or not isclose(n_in, n_out, rel_tol=1e-8):
        return None
    expected = e_in * (n_out / n_in)
    if isclose(expected, e_out, rel_tol=1e-6, abs_tol=.001):
        return None
    origin, target = first.get("reference_option_hints", {}), last.get("reference_option_hints", {})
    changed = []
    for field, label in (("section", "section"), ("producer_free_water", "producer BLKOPFREWAT"),
                         ("producer_phase_options", "producer phase settings"),
                         ("producer_property_options", "producer property settings")):
        a, b = origin.get(field), target.get(field)
        if a != b:
            changed.append(f"{label} {a or 'undeclared'} to {b or 'undeclared'}")
    details = ", ".join(changed) or f"option keys {keys[0]} and {keys[1]}"
    return _finding("Cross-option chemical exergy difference", "block", block["id"],
                    "outlet_minus_scaled_inlet_chemical_exergy_kw", e_out - expected, "kW",
                    f"{first['id']} to {last['id']} has matching composition and throughput under {method}. "
                    f"Saved settings differ ({details}). The value compares chemical exergy at outlet throughput.",
                    "Check separate option-matched Aspen ambient flashes for both streams, including free-water and allowed phases.")


def build_quality(report):
    """Return review records with signed values on the basis of the supplied analysis.

    Keep five independent-model differences with the largest absolute magnitudes.
    Remaining differences are available per stream. Balance flags use the
    tolerances of the original control-volume calculation.
    """
    rows = []
    streams = report["streams"].values()
    comparisons = sorted((s for s in streams if s.get("independent_difference_kw")),
                         key=lambda s: (-abs(s["independent_difference_kw"]), s["id"]))
    for s in comparisons[:5]:
        rows.append(_finding("Physical exergy model difference", "stream", s["id"], "independent_difference_kw",
                             s["independent_difference_kw"], "kW",
                             f"Physical exergy using anchor {s['reference_anchor_id']} minus the independent estimate. "
                             f"Ranked by absolute difference across {len(comparisons)} nonzero comparisons.",
                             "Check the stream and anchor with Aspen EXERGYFL or a common-composition ambient flash."))
    cross_option = [item for block in report["blocks"]
                    if (item := _cross_option_comparison(block, report["streams"])) is not None]
    rows.extend(sorted(cross_option, key=lambda row: (-abs(row["value"]), row["id"]))[:5])
    for t in report["transfers"]:
        key = f"{t['block']}/{t['id']}"
        if t.get("rate_kw") is None:
            rows.append(_finding("Missing transfer rate", "transfer", key, "rate_kw", None, "kW",
                                 f"The {t['kind']} transfer has no saved rate.",
                                 "Supply the signed rate from the utility or shaft-work record."))
        elif t["kind"] == "heat" and t["rate_kw"] != 0 and t.get("temperature_k") is None:
            rows.append(_finding("Missing heat boundary temperature", "transfer", key, "temperature_k", None, "K",
                                 f"The signed heat rate is {t['rate_kw']:+.6g} kW.",
                                 "Set temperature_k in the transfer settings from the utility contact conditions."))
    groups = list(report.get("groups", {}).values())
    for scope, boundaries in (("block", report["blocks"]), ("group", groups), ("plant", [report["plant"]])):
        for boundary in boundaries:
            rows.extend(_boundary_findings(scope, boundary, groups))
    for stream in streams:
        rows.extend(_stream_findings(stream))
    equipment = {b["id"] for b in report["blocks"]}
    for warning in report.get("source_warnings", []):
        text = str(warning)
        origin, _, detail = text.partition(": ")
        if origin in equipment:
            scope, item = "block", origin
        elif origin in report["streams"]:
            scope, item = "stream", origin
        else:
            scope, item, detail = "source", report["source"].get("source_filename", "case"), text
        rows.append(_finding("Saved Aspen warning", scope, item, "Aspen warning", None, "", detail,
                             "Review the original warning and saved convergence or phase-stability results."))
    return rows
