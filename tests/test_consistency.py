"""Scientific reporting checks: disjoint boundaries and unmodified native data."""
from copy import deepcopy
import unittest

from aspen_energy_exergy.consistency import build_consistency_audit


def fixture():
    def stream(name, chemical, options, free_water):
        return {
            "id": name, "molar_kmol_h": 3600.0, "component_flows_kmol_h": {"H2O": 3600.0},
            "enthalpy_kw": -10000.0, "entropy_kw_k": -3.0,
            "physical_exergy_kw": 10.0, "chemical_exergy_kw": chemical,
            "total_exergy_kw": chemical + 10.0, "property_method": "RKS-BM",
            "reference_options_key": options,
            "reference_option_hints": {"producer_free_water": free_water,
                                       "section_options": {"GFREE-WATER": "STEAMNBS"}},
        }
    blocks = [
        {"id": "A", "exergy_destruction_kw": -2.0, "energy_residual_kw": 0.0},
        {"id": "B", "exergy_destruction_kw": 5.0, "energy_residual_kw": 0.0},
        {"id": "C", "exergy_destruction_kw": 1.0, "energy_residual_kw": 0.0},
    ]
    return {
        "streams": {"WAT": stream("WAT", 727.0, "free-water", "YES"),
                    "WAT-C": stream("WAT-C", 184.0, "eos-water", None)},
        "blocks": blocks,
        "groups": {
            "guard_beds": {"members": ["A", "B"], "exergy_destruction_kw": 3.0,
                           "energy_residual_kw": 0.0},
            "all_membranes": {"members": ["A", "B"], "exergy_destruction_kw": 3.0},
        },
        "plant": {"exergy_destruction_kw": 5.0},
        "transfers": [
            {"link_id": "Q", "block": "B", "rate_kw": -10.0, "exergy_kw": -2.0},
            {"link_id": "Q", "block": "C", "rate_kw": 10.0, "exergy_kw": 1.0},
        ],
        "dryer_design": {"regeneration_heat_kw": 120.0, "regeneration_heat_exergy_kw": 47.4500659029,
                         "required_mass_per_bed_kg": 5543.34411489, "installed_mass_per_bed_kg": 5540.0},
    }


class ConsistencyAuditTests(unittest.TestCase):
    def test_partition_counts_each_equipment_once_and_retains_cross_link(self):
        report = fixture()
        result = build_consistency_audit(report)
        part = result["package_reporting"]
        self.assertEqual(part["equipment_owner"], {"A": "guard_beds", "B": "guard_beds", "C": "C"})
        self.assertEqual(part["totals"]["owner_count"], 2)
        self.assertEqual(part["totals"]["owner_destruction_sum_kw"], 4.0)
        self.assertEqual(part["totals"]["cross_owner_link_destruction_kw"], 1.0)
        self.assertEqual(part["totals"]["reconstructed_plant_destruction_kw"], 5.0)
        self.assertTrue(part["totals"]["reconciles"])
        self.assertEqual(part["excluded_overlapping_views"], ["all_membranes"])

    def test_overlap_is_rejected_instead_of_double_counted(self):
        report = fixture()
        report["groups"]["psa"] = {"members": ["B"], "exergy_destruction_kw": 5.0}
        with self.assertRaisesRegex(ValueError, "overlapping"):
            build_consistency_audit(report)

    def test_missing_link_or_boundary_exergy_is_never_zero_filled(self):
        report = fixture()
        report["transfers"][1]["exergy_kw"] = None
        totals = build_consistency_audit(report)["package_reporting"]["totals"]
        self.assertIsNone(totals["cross_owner_link_destruction_kw"])
        self.assertIsNone(totals["reconstructed_plant_destruction_kw"])
        self.assertIsNone(totals["reconciles"])
        report = fixture()
        report["blocks"][2]["exergy_destruction_kw"] = None
        self.assertIsNone(build_consistency_audit(report)["package_reporting"]["totals"]["owner_destruction_sum_kw"])

    def test_water_option_change_is_not_silently_standardized(self):
        report = fixture()
        water = build_consistency_audit(report)["water_reference"]
        self.assertEqual(len(water["families"]), 2)
        transition = water["WAT_to_WAT_C"]
        self.assertFalse(transition["same_reference_options"])
        self.assertEqual(transition["inlet"]["effective_water_method"], "STEAMNBS")
        self.assertEqual(transition["outlet"]["effective_water_method"], "RKS-BM")
        self.assertEqual(transition["inlet_minus_outlet_chemical_kw"], 543.0)
        self.assertEqual(transition["status"], "different_water_property_bases")
        for row in report["streams"].values():
            row.pop("reference_options_key")
        water = build_consistency_audit(report)["water_reference"]
        self.assertIsNone(water["WAT_to_WAT_C"]["same_reference_options"])
        self.assertTrue(all(g["same_option_reference_consistent"] is None for g in water["families"]))

    def test_negative_values_and_resources_preserved_without_mutation(self):
        report = fixture()
        before = deepcopy(report)
        result = build_consistency_audit(report, previous=deepcopy(report))
        self.assertEqual(report, before)
        negative = result["package_reporting"]["negative_block_diagnostics"][0]
        self.assertEqual(negative["raw_exergy_destruction_kw"], -2.0)
        self.assertEqual(negative["owner_exergy_destruction_kw"], 3.0)
        dryer = result["resources"]["dryer"]
        self.assertEqual(dryer["regeneration_heat_exergy_kw"], 47.4500659029)
        self.assertEqual(dryer["minimum_integer_sieve_mass_per_bed_kg"], 5544)
        self.assertEqual(dryer["added_to_plant_exergy_destruction_kw"], 0.0)
        self.assertEqual(dryer["added_to_plant_electricity_kw"], 0.0)
        self.assertEqual(result["resources"]["psa"]["added_to_native_balance_kw"], 0.0)
        self.assertFalse(result["unmodeled_resources_added"])

    def test_secondary_missing_split_is_explicit_not_a_primary_gap(self):
        report = fixture()
        saved = deepcopy(report["streams"])
        saved["WAT"]["physical_exergy_kw"] = None
        saved["WAT"]["chemical_exergy_kw"] = None
        report["archive_state_analysis"] = {
            "streams": saved,
            "reference_availability": {"counts": {"unavailable_changed_recipe": 1},
                                       "streams": [{"id": "WAT", "status": "unavailable_changed_recipe"}]},
        }
        result = build_consistency_audit(report)
        scope = result["state_scope"]
        self.assertEqual(scope["primary"]["physical_chemical_splits_available"], 2)
        self.assertEqual(scope["saved_archive"]["physical_chemical_splits_available"], 1)
        self.assertEqual(scope["saved_archive"]["total_exergy_available"], 2)
        self.assertEqual(result["export_decision"]["new_state_matched_reference_needed_for_secondary_ids"], ["WAT"])
        self.assertFalse(result["export_decision"]["source_recovery_review_complete"])

    def test_native_recovery_resolves_current_request_without_backfilling_history(self):
        report=fixture()
        report["native_rerun"]={"completed":True}
        saved=deepcopy(report["streams"])
        saved["WAT"]["physical_exergy_kw"]=None
        saved["WAT"]["chemical_exergy_kw"]=None
        report["archive_state_analysis"]={"streams":saved,"reference_availability":{
            "streams":[{"id":"WAT","status":"unavailable_changed_recipe"}]}}
        result=build_consistency_audit(report)
        self.assertFalse(result["export_decision"]["new_user_export_required"])
        self.assertTrue(result["export_decision"]["historical_gaps_resolved_in_current_snapshot"])
        self.assertEqual(result["export_decision"]["new_state_matched_reference_needed_for_secondary_ids"],[])
        self.assertEqual(result["export_decision"]["historical_reference_gap_ids"],["WAT"])
        self.assertIsNone(report["archive_state_analysis"]["streams"]["WAT"]["physical_exergy_kw"])
        report["streams"]["WAT"]["physical_exergy_kw"]=None
        result=build_consistency_audit(report)
        self.assertTrue(result["export_decision"]["new_user_export_required"])
        self.assertFalse(result["export_decision"]["historical_gaps_resolved_in_current_snapshot"])


if __name__ == "__main__":
    unittest.main()
