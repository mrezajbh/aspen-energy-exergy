"""Auditable reporting boundaries and remaining Aspen data requirements.

This module does not calculate replacement thermodynamic properties, modify an
Aspen state, or add unmodeled equipment resources. It organizes the calculated
native-engine results and exposes consistency limitations separately.
"""
from collections import defaultdict
from math import ceil, fsum, isclose, isfinite


PRIMARY_PACKAGES = (
    "guard_beds", "biogas_membranes", "psa", "co_membrane", "dryer",
    "ft_heat_recovery", "cooler_water_section",
)
OVERLAPPING_VIEWS = ("all_membranes", "selected_packages")


def _number(value):
    return value if type(value) in (float, int) and isfinite(value) else None


def _sum(values):
    values = [_number(v) for v in values]
    return None if None in values else fsum(values)


def _difference(first, second):
    first, second = _number(first), _number(second)
    return None if first is None or second is None else first - second


def _coverage(streams):
    return {
        "streams": len(streams),
        "enthalpy_available": sum(_number(s.get("enthalpy_kw")) is not None for s in streams.values()),
        "entropy_available": sum(_number(s.get("entropy_kw_k")) is not None for s in streams.values()),
        "total_exergy_available": sum(_number(s.get("total_exergy_kw")) is not None for s in streams.values()),
        "physical_chemical_splits_available": sum(
            _number(s.get("physical_exergy_kw")) is not None and
            _number(s.get("chemical_exergy_kw")) is not None for s in streams.values()),
    }


def _effective_water_method(stream):
    hints = stream.get("reference_option_hints") or {}
    options = hints.get("section_options") or {}
    if str(hints.get("producer_free_water", "")).upper() == "YES":
        return options.get("GFREE-WATER", "free-water method not recorded")
    return stream.get("property_method")


def _water_audit(streams):
    rows = []
    grouped = defaultdict(list)
    for sid, stream in streams.items():
        flows = stream.get("component_flows_kmol_h") or {}
        positive = {name for name, value in flows.items() if _number(value) is not None and value > 0}
        n = _number(stream.get("molar_kmol_h"))
        if positive != {"H2O"} or n is None or n <= 0:
            continue
        chemical = _number(stream.get("chemical_exergy_kw"))
        key = stream.get("reference_options_key")
        method = stream.get("property_method")
        hints = stream.get("reference_option_hints") or {}
        row = {
            "id": sid, "molar_kmol_h": n, "property_method": method,
            "effective_water_method": _effective_water_method(stream),
            "reference_options_key": key,
            "producer_free_water": hints.get("producer_free_water"),
            "option_metadata_available": bool(key and hints),
            "chemical_exergy_j_mol": None if chemical is None else chemical * 3600.0 / n,
            "physical_exergy_kw": _number(stream.get("physical_exergy_kw")),
            "chemical_exergy_kw": chemical,
            "temperature_k": _number(stream.get("temperature_k")),
            "pressure_pa": _number(stream.get("pressure_pa")),
            "enthalpy_kw": _number(stream.get("enthalpy_kw")),
            "entropy_kw_k": _number(stream.get("entropy_kw_k")),
        }
        rows.append(row)
        # Missing provenance is not permission to treat two option sets as equal.
        family = (method, key) if key and hints else (method, "unverified:" + sid)
        grouped[family].append(row)
    families = []
    for (method, key), members in grouped.items():
        values = [r["chemical_exergy_j_mol"] for r in members
                  if r["chemical_exergy_j_mol"] is not None]
        span = max(values) - min(values) if values else None
        verified = not key.startswith("unverified:")
        families.append({
            "property_method": method, "reference_options_key": key if verified else None,
            "metadata_verified": verified, "streams": [r["id"] for r in members],
            "available_chemical_values": len(values),
            "chemical_exergy_min_j_mol": min(values) if values else None,
            "chemical_exergy_max_j_mol": max(values) if values else None,
            "within_family_range_j_mol": span,
            "same_option_reference_consistent": None if not verified or len(values) < 2 else span <= 1e-4,
        })
    lookup = {r["id"]: r for r in rows}
    a, b = lookup.get("WAT"), lookup.get("WAT-C")
    transition = {"available": bool(a and b), "states_modified": False,
                  "requires_consistent_state_review": False}
    if a and b:
        same_flow = isclose(a["molar_kmol_h"], b["molar_kmol_h"], rel_tol=1e-9, abs_tol=1e-9)
        verified = a["option_metadata_available"] and b["option_metadata_available"]
        same_options = (a["property_method"] == b["property_method"] and
                        a["reference_options_key"] == b["reference_options_key"]) if verified else None
        chemical_difference = _difference(a["chemical_exergy_kw"], b["chemical_exergy_kw"])
        transition.update({
            "inlet": a, "outlet": b, "same_pure_water_throughput": same_flow,
            "same_reference_options": same_options,
            "inlet_minus_outlet_chemical_kw": chemical_difference,
            "requires_consistent_state_review": chemical_difference is None or (same_flow and abs(chemical_difference) > 1e-6),
            "inlet_minus_outlet_physical_kw": _difference(a["physical_exergy_kw"], b["physical_exergy_kw"]),
            "native_enthalpy_drop_kw": _difference(a["enthalpy_kw"], b["enthalpy_kw"]),
            "status": "different_water_property_bases" if same_options is False else "reference_consistency_requires_review",
            "interpretation": "A basis/interface diagnostic, not a chemical reaction or verified physical irreversibility. Raw Aspen H/S and EXERGYFL are retained.",
            "needs": "Operating and ambient states recalculated with a consistent water-property/free-water convention; changing only the assigned chemical reference is not a complete thermodynamic correction.",
        })
    return {
        "streams": rows, "families": families, "WAT_to_WAT_C": transition,
        "within_family_tolerance_j_mol": 1e-4,
        "classification_basis": "Strict property-method plus saved reference-option identity. Effective free-water method is shown separately. Different option families are not called same-reference conflicts.",
    }


def _package_audit(report, previous):
    blocks = {b["id"]: b for b in report["blocks"]}
    if len(blocks) != len(report["blocks"]):
        raise ValueError("Duplicate equipment identifier in consistency audit")
    groups = report.get("groups", {})
    old_blocks = {b["id"]: b for b in (previous or {}).get("blocks", [])}
    old_groups = (previous or {}).get("groups", {})
    partition, owner = [], {}
    for name in PRIMARY_PACKAGES:
        if name not in groups:
            continue
        group = groups[name]
        members = group.get("members", [])
        if not members or len(members) != len(set(members)):
            raise ValueError(f"Package {name} has empty or duplicate members")
        for bid in members:
            if bid not in blocks or bid in owner:
                raise ValueError(f"Unknown or overlapping primary package member: {bid}")
            owner[bid] = name
        old = old_groups.get(name, {})
        prior = old.get("exergy_destruction_kw") if set(old.get("members", [])) == set(members) else None
        partition.append({
            "id": name, "kind": "package", "members": list(members),
            "exergy_destruction_kw": _number(group.get("exergy_destruction_kw")),
            "energy_residual_kw": _number(group.get("energy_residual_kw")),
            "work_input_kw": _number(group.get("work_input_kw")),
            "heat_input_kw": _number(group.get("heat_input_kw")),
            "heat_rejected_kw": _number(group.get("heat_rejected_kw")),
            "previous_exergy_destruction_kw": _number(prior),
            "change_kw": _difference(group.get("exergy_destruction_kw"), prior),
            "qualification": "Conditional water-property interface diagnostic" if name == "cooler_water_section"
            else "Modeled package boundary; omitted consumables/regeneration remain separately identified",
        })
    for bid, block in blocks.items():
        if bid in owner:
            continue
        owner[bid] = bid
        prior = old_blocks.get(bid, {}).get("exergy_destruction_kw")
        partition.append({
            "id": bid, "kind": "equipment", "members": [bid],
            "exergy_destruction_kw": _number(block.get("exergy_destruction_kw")),
            "energy_residual_kw": _number(block.get("energy_residual_kw")),
            "previous_exergy_destruction_kw": _number(prior),
            "change_kw": _difference(block.get("exergy_destruction_kw"), prior),
            "qualification": "Imported-state equipment ledger; native-equation coverage is reported separately",
        })
    legs = defaultdict(list)
    for transfer in report.get("transfers", []):
        if transfer.get("link_id") is not None:
            legs[transfer["link_id"]].append(transfer)
    cross_links = []
    for name, pair in legs.items():
        if len(pair) != 2 or any(t.get("block") not in owner for t in pair):
            raise ValueError(f"Incomplete transfer-link inventory: {name}")
        owners = [owner[t["block"]] for t in pair]
        if owners[0] == owners[1]:
            continue  # Already included in that package's external boundary.
        signed = _sum(t.get("exergy_kw") for t in pair)
        cross_links.append({
            "id": name, "owners": owners, "blocks": [t["block"] for t in pair],
            "exergy_destruction_kw": None if signed is None else -signed,
            "energy_cancellation_residual_kw": _sum(t.get("rate_kw") for t in pair),
            "basis": "Negative signed sum of two calculated transfer legs; no extra equipment duty",
        })
    unit_sum = _sum(r["exergy_destruction_kw"] for r in partition)
    link_sum = _sum(r["exergy_destruction_kw"] for r in cross_links)
    total = _sum((unit_sum, link_sum))
    plant = _number(report["plant"].get("exergy_destruction_kw"))
    residual = _difference(total, plant)
    diagnostic = []
    lookup = {r["id"]: r for r in partition}
    for bid, block in blocks.items():
        value = _number(block.get("exergy_destruction_kw"))
        if value is not None and value < -1e-3:
            own = owner[bid]
            diagnostic.append({"id": bid, "raw_exergy_destruction_kw": value,
                               "reporting_owner": own, "owner_kind": lookup[own]["kind"],
                               "owner_exergy_destruction_kw": lookup[own]["exergy_destruction_kw"],
                               "negative_value_preserved": True})
    return {
        "partition": partition, "equipment_owner": owner,
        "negative_block_diagnostics": diagnostic,
        "cross_owner_transfer_links": cross_links,
        "excluded_overlapping_views": [k for k in OVERLAPPING_VIEWS if k in groups],
        "totals": {
            "equipment_count": len(blocks), "owner_count": len(partition),
            "package_count": sum(r["kind"] == "package" for r in partition),
            "singleton_count": sum(r["kind"] == "equipment" for r in partition),
            "equipment_coverage_count": len(owner), "overlap_count": 0,
            "owner_destruction_sum_kw": unit_sum, "cross_owner_link_destruction_kw": link_sum,
            "reconstructed_plant_destruction_kw": total, "plant_destruction_kw": plant,
            "partition_minus_plant_kw": residual,
            "reconciles": None if residual is None else isclose(total, plant, rel_tol=1e-9, abs_tol=1e-6),
        },
        "interpretation": "Each equipment block has exactly one reporting owner. Nested comparison groups are excluded from totals. Positive package destruction does not validate conceptual equipment physics or remove property-reference artifacts.",
    }


def _resource_audit(report):
    supplied = report.get("dryer_design", {})
    fields = ("water_removed_kg_h", "bed_count", "installed_mass_per_bed_kg", "required_mass_per_bed_kg",
              "bed_mass_shortfall_kg", "working_capacity_kg_per_kg", "adsorption_h", "regeneration_h",
              "purge_co2_kmol_h", "regeneration_heat_kw", "regeneration_temperature_c", "regeneration_heat_exergy_kw")
    dryer = {key: supplied.get(key) for key in fields}
    required = _number(supplied.get("required_mass_per_bed_kg"))
    dryer.update({
        "minimum_integer_sieve_mass_per_bed_kg": None if required is None else ceil(required),
        "heat_is_electricity": False, "included_in_native_balance": False,
        "added_to_plant_exergy_destruction_kw": 0.0, "added_to_plant_electricity_kw": 0.0,
        "resource_balance_complete": False,
        "interpretation": "Supplied regeneration heat exergy is an input allowance, not additional destruction. Continuous alternating regeneration would require the stated mean duty continuously; peak/hold/cool timing is not established.",
        "missing": ["Regeneration purge inlet/outlet T, P, flow and water content",
                    "Whether heat is peak, regeneration-period mean or whole-cycle mean",
                    "Utility origin, recycle destination and overlap with existing Aspen utilities",
                    "Cycle heat/entropy closure including regeneration cooling and media state"],
    })
    psa = report.get("groups", {}).get("psa", {})
    return {
        "dryer": dryer,
        "psa": {
            "modeled_package_work_input_kw": _number(psa.get("work_input_kw")),
            "report_design_electricity_kw": 900.0,
            "report_design_cooling_kw": 900.0,
            "added_to_native_balance_kw": 0.0,
            "basis": "BtJ Final Report section 3.3.5 assigns a separate 0.9 MW PSA design allowance; the user authorized K-301, LOWP and K-802 as package resources.",
            "incremental_to_existing_compressors": None,
            "interpretation": "Keep the report allowance separate; its overlap with existing package compression is not resolved. Do not add both as if their independence were demonstrated.",
        },
        "guard_beds": {
            "fresh_spent_sorbents_modeled": False,
            "added_to_native_balance_kw": 0.0,
            "known_basis": "Report specifies ZnO + H2S -> ZnS + H2O and activated-carbon loading of 10-20 wt%.",
            "missing": ["Fresh/spent sorbent rates, composition and thermodynamic states",
                        "Replacement/regeneration basis and associated heat/work flows"],
            "interpretation": "Existing gas-side package is retained; a complete resource exergy cannot be inferred from a conceptual impurity outlet alone.",
        },
    }


def build_consistency_audit(report, previous=None):
    """Return a JSON-ready audit without mutating ``report`` or ``previous``.

    The Biogas-to-Jet reporting groups are optional for smaller test cases, but
    their primary members must be disjoint whenever present. Existing native
    engine outputs are used directly. The final export decision can be refined
    by the calling source-recovery workflow; this function does not claim that
    an unsearched source file is absent.
    """
    streams = report["streams"]
    archive = report.get("archive_state_analysis", {})
    saved_streams = archive.get("streams", {})
    reference_rows = archive.get("reference_availability", {}).get("streams", [])
    unavailable = [r["id"] for r in reference_rows if r.get("status") == "unavailable_changed_recipe"]
    primary = _coverage(streams)
    saved = _coverage(saved_streams)
    water = _water_audit(streams)
    fresh = report.get("native_rerun", {}).get("completed") is True
    primary_unavailable = [sid for sid, stream in streams.items()
                           if _number(stream.get("physical_exergy_kw")) is None or
                           _number(stream.get("chemical_exergy_kw")) is None]
    recovered_ids = [sid for sid in unavailable if sid in streams and sid not in primary_unavailable] if fresh else []
    return {
        "schema": "aspen-consistency-audit/v1",
        "states_modified": False, "unmodeled_resources_added": False,
        "state_scope": {
            "primary": primary, "saved_archive": saved,
            "saved_archive_split_counts": archive.get("reference_availability", {}).get("counts", {}),
            "saved_archive_unmatched_split_ids": unavailable,
            "saved_archive_is_historical": fresh,
            "primary_unmatched_split_ids": primary_unavailable,
            "historical_gap_ids_available_in_current_snapshot": recovered_ids,
            "primary_basis": report.get("comparison_basis", {}).get("primary"),
            "interpretation": ("The fresh native Aspen snapshot is primary. Former saved-state gaps are retained as historical provenance; current native values are not backfilled into that older state." if fresh else
                               "Missing secondary splits do not mean primary stream energy or total exergy is missing. State sources are not silently mixed."),
        },
        "water_reference": water,
        "package_reporting": _package_audit(report, previous),
        "resources": _resource_audit(report),
        "export_decision": {
            "source_recovery_review_complete": fresh,
            "new_export_required_for_existing_primary_ledger": False,
            "new_state_matched_reference_needed_for_secondary_ids": [] if fresh else unavailable,
            "historical_reference_gap_ids": unavailable,
            "current_primary_reference_gap_ids": primary_unavailable,
            "historical_gaps_resolved_in_current_snapshot": len(recovered_ids) == len(unavailable) if fresh else False,
            "new_user_export_required": bool(primary_unavailable) if fresh else None,
            "secondary_reference_alternatives": [
                "Native EXERGYFL matching the saved operating composition/options and substream scope",
                "Ambient H0 and S0 for exactly the same composition, substreams and property/free-water options",
            ],
            "water_interface_needs_consistent_operating_and_ambient_states": water["WAT_to_WAT_C"]["requires_consistent_state_review"],
            "conceptual_standalone_resources_need_engineering_data_not_just_an_export": True,
            "interpretation": ("A fresh native Aspen run recovered the current snapshot from the existing project. No additional user export is required for the primary ledger. Historical saved-state reference gaps remain documented in their original snapshot. Water-interface physical validity and omitted engineering resources remain separate issues." if fresh and not primary_unavailable else
                               "Existing primary data support imported-state/package reporting. Complete secondary splits and validated individual water-interface balances require matched reference evidence; the caller must record whether that evidence was recovered from existing files or remains to be exported."),
        },
    }
