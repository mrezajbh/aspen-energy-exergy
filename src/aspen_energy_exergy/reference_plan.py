"""Plan same-composition Aspen equilibrium flashes at the reference environment."""

from copy import deepcopy
from math import isclose, log

from ._exchange import case_hashes, digest
from .properties import P0_PA, REFERENCE_PLAN_SCHEMA, T0_K, calculate_streams, equivalent_reference_streams


def _at_environment(stream):
    t, p = stream.get("temperature_k"), stream.get("pressure_pa")
    return (t is not None and p is not None
            and isclose(t, T0_K, abs_tol=1e-6, rel_tol=0.)
            and isclose(p, P0_PA, abs_tol=1e-4, rel_tol=0.))


def _representative_order(stream):
    t, p = stream.get("temperature_k"), stream.get("pressure_pa")
    distance = abs(log(t / T0_K)) + abs(log(p / P0_PA)) if t and p and min(t, p) > 0 else float("inf")
    return (not _at_environment(stream), distance, stream["id"])


def _groups(streams):
    """Fix each representative before admitting peers, including near-equal mixtures."""
    groups = []
    for stream in sorted(streams, key=_representative_order):
        for members in groups:
            if equivalent_reference_streams(members[0], stream):
                members.append(stream)
                break
        else:
            groups.append([stream])
    return sorted(groups, key=lambda members: (-len(members), members[0]["id"]))


def _missing_input(stream):
    if not stream.get("property_method"):
        return "Supply the operating stream's Aspen property method."
    if not stream.get("molar_kmol_h") or not any((stream.get("component_flows_kmol_h") or {}).values()):
        return "Supply the positive molar throughput and component rates for the operating stream."
    return None


def _candidate_ids(representative, candidates):
    return sorted(s["id"] for s in candidates
                  if s["id"] == representative["id"] or equivalent_reference_streams(representative, s))


def _record(members, component_ids, ambient, existing):
    rep = members[0]
    ids = [rep["id"], *sorted(s["id"] for s in members[1:])]
    subs = rep.get("substreams") or {}
    subflows = {name: deepcopy(sub.get("component_flows_kmol_h")) for name, sub in sorted(subs.items())}
    phase_sets = sorted({tuple(sorted(s.get("active_phases") or [])) for s in members})
    notes = []
    if not rep.get("property_options_id"):
        notes.append("Verify the complete Aspen property and phase settings. The saved method and any local option declarations provide partial evidence.")
    if len(phase_sets) > 1:
        notes.append("These streams occupy different operating phases. Verify equivalent allowed phases at the ambient flash.")
    if len(subs) > 1:
        notes.append("Preserve the component allocation across each active substream and report whole-material H and S.")
    if not subs:
        notes.append("Check the substream inventory for this normalized export in Aspen.")
    member_sources = {}
    for s in members:
        hints = s.get("reference_option_hints") or {}
        member_sources[s["id"]] = {
            "producer_block": s.get("from_block"), "section": hints.get("section"),
            "exchanger_side": hints.get("exchanger_side"), "active_phases": list(s.get("active_phases") or []),
        }
    return {
        "id": rep["id"], "operating_state_sha256": digest(rep),
        "member_ids": ids, "member_count": len(ids), "applies_to": ids[1:],
        "member_sources": member_sources,
        "property_method": rep["property_method"], "property_options_id": rep.get("property_options_id"),
        "reference_options_key": rep.get("reference_options_key"),
        "reference_option_hints": deepcopy(rep.get("reference_option_hints")),
        "operating_temperature_k": rep.get("temperature_k"), "operating_pressure_pa": rep.get("pressure_pa"),
        "operating_phase_sets": [list(phases) for phases in phase_sets],
        "substream_component_flows_kmol_h": subflows,
        "saved_ambient_candidate_ids": _candidate_ids(rep, ambient),
        "candidate_existing_flash_ids": _candidate_ids(rep, existing),
        "notes": notes,
        "component_ids": list(component_ids),
        "component_flows_kmol_h": [rep["component_flows_kmol_h"].get(cid, 0.) for cid in component_ids],
        "T_K": T0_K, "P_Pa": P0_PA, "phase_convention": "equilibrium", "quantity_basis": "flow_SI",
        "same_aspen_property_options": False, "same_aspen_property_options_for_all": False,
        "aspen_reference_id": None, "source_location": None, "H0_kW": None, "S0_kW_K": None,
    }


def build_reference_plan(model, *, references=None, native_exergy=None):
    """Return a fillable reference document and an audit of the pending streams.

    Supply a loaded case, an optional previously returned reference document,
    and optional validated native EXERGYFL records. This pass uses only saved
    results. Member lists identify candidate reuse subject to operator checks.
    """
    streams = model["streams"]
    rows = calculate_streams(model, references=references, native_exergy=native_exergy, estimate_splits=False)
    covered, unplannable, pending, ambient, existing = {}, [], [], [], []
    for stream in streams:
        sid = stream["id"]
        result = rows[sid]
        if (stream.get("molar_kmol_h") and _at_environment(stream)
                and all(stream.get(k) is not None for k in ("enthalpy_kw", "entropy_kw_k"))):
            ambient.append(stream)
        if result["split_method"] == "aspen_ambient_flash":
            existing.append(stream)
        if result["split_quality"] != "unavailable":
            covered[sid] = result["split_quality"]
            continue
        reason = _missing_input(stream)
        if reason:
            unplannable.append({"id": sid, "reason": reason})
        else:
            pending.append(stream)

    component_ids = [c["id"] for c in model["components"]]
    records = [_record(members, component_ids, ambient, existing) for members in _groups(pending)]
    return {
        "plan_schema": REFERENCE_PLAN_SCHEMA, "case_sha256": case_hashes(model)[0],
        "source_filename": model["metadata"].get("source_filename"),
        "environment": {"temperature_k": T0_K, "pressure_pa": P0_PA},
        "summary": {
            "total_streams": len(streams), "already_covered_streams": len(covered),
            "pending_streams": len(streams) - len(covered), "plannable_streams": len(pending),
            "unplannable_streams": len(unplannable), "flash_groups": len(records),
            "saved_ambient_candidates": len(ambient),
        },
        "instructions": [
            "Rows with more potential users appear first. Representatives are selected by proximity to the reference temperature and pressure.",
            "Flash each representative composition at T_K and P_Pa with the same Aspen property method, component data, property options, allowed equilibrium phases and substream allocation.",
            "Enter whole-stream H0_kW and S0_kW_K for the component flows listed, plus the Aspen reference ID and source location. Verify the source settings before marking same_aspen_property_options true.",
            "Before reusing the flash, retain only verified peer IDs in applies_to and set same_aspen_property_options_for_all true. For a representative-only result use an empty applies_to list.",
            "A saved stream at the ambient temperature and pressure is a candidate only. Verify the equilibrium and reference settings in Aspen before returning its values.",
        ],
        "already_covered": dict(sorted(covered.items())),
        "unplannable": sorted(unplannable, key=lambda item: item["id"]),
        "saved_ambient_candidates": [
            {k: s.get(k) for k in ("id", "property_method", "temperature_k", "pressure_pa", "molar_kmol_h",
                                   "enthalpy_kw", "entropy_kw_k")}
            for s in sorted(ambient, key=lambda s: s["id"])
        ],
        "streams": records,
    }
