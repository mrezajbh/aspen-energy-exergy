import copy
from pathlib import Path
import unittest

from aspen_energy_exergy.api import load_case
from aspen_energy_exergy.properties import calculate_streams
from aspen_energy_exergy.reference_plan import build_reference_plan


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "outputs/energy_exergy_aspen_20261111/source/Biogas to Jetfuel - FT.apwz"


def heater():
    model = load_case(ROOT / "examples/heater.json")
    for s in model["streams"]:
        s["substreams"] = {"MIXED": {"component_flows_kmol_h": dict(s["component_flows_kmol_h"]),
                                    "molar_kmol_h": s["molar_kmol_h"]}}
    return model


class ReferencePlanTest(unittest.TestCase):
    def test_ambient_candidate_is_unconfirmed_until_the_operator_returns_it(self):
        model = heater()
        before = copy.deepcopy(model)
        plan = build_reference_plan(model)
        self.assertEqual(before, model)
        self.assertEqual((plan["summary"]["flash_groups"], plan["summary"]["saved_ambient_candidates"]), (1, 1))
        row, = plan["streams"]
        self.assertEqual(row["id"], "FEED")
        self.assertEqual(row["saved_ambient_candidate_ids"], ["FEED"])
        self.assertEqual(row["applies_to"], ["PRODUCT"])
        self.assertEqual(row["substream_component_flows_kmol_h"], {"MIXED": {"H2": 1}})
        self.assertIsNone(row["H0_kW"])
        self.assertFalse(row["same_aspen_property_options"])
        self.assertFalse(row["same_aspen_property_options_for_all"])
        raw = calculate_streams(model, estimate_splits=False, references=plan)
        self.assertTrue(all(r["physical_exergy_kw"] is None for r in raw.values()))

        row.update(aspen_reference_id="Reference 1", source_location="Synthetic heater equilibrium case",
                   H0_kW=0., S0_kW_K=0., same_aspen_property_options=True,
                   same_aspen_property_options_for_all=True)
        resolved = calculate_streams(model, estimate_splits=False, references=plan)
        self.assertEqual(resolved["FEED"]["physical_exergy_kw"], 0.)
        self.assertAlmostEqual(resolved["PRODUCT"]["physical_exergy_kw"],
                               0.8 - 298.15 * 0.0023138568250747125)
        self.assertEqual(build_reference_plan(model, references=plan)["streams"], [])

    def test_prior_direct_flash_is_suggested_for_an_unconfirmed_peer(self):
        model = heater()
        returned = build_reference_plan(model)
        row, = returned["streams"]
        row.update(aspen_reference_id="Reference 1", source_location="Synthetic fixture",
                   H0_kW=0., S0_kW_K=0., same_aspen_property_options=True, applies_to=[])
        updated = build_reference_plan(model, references=returned)
        self.assertEqual(updated["already_covered"], {"FEED": "matched"})
        self.assertEqual(updated["streams"][0]["candidate_existing_flash_ids"], ["FEED"])
        self.assertEqual(updated["streams"][0]["member_ids"], ["PRODUCT"])

    def test_different_options_and_unknown_substream_allocations_stay_individual(self):
        for key, alternative in (("property_method", "OTHER"), ("property_options_id", "section B"),
                                 ("reference_options_key", "free-water:yes"),
                                 ("substreams", {}), ("substreams", {"CISOLID": {"component_flows_kmol_h": {"H2": 1}}})):
            with self.subTest(key=key, alternative=alternative):
                model = heater()
                model["streams"][1][key] = alternative
                rows = build_reference_plan(model)["streams"]
                self.assertEqual(len(rows), 2)
                self.assertTrue(all(r["applies_to"] == [] for r in rows))

    def test_multiple_substreams_require_matching_allocation_not_just_matching_total(self):
        model = heater()
        for s, x in zip(model["streams"], (0.3, 0.4)):
            s["substreams"] = {name: {"component_flows_kmol_h": {"H2": part}, "molar_kmol_h": part}
                               for name, part in (("MIXED", x), ("CISOLID", 1 - x))}
        self.assertEqual(build_reference_plan(model)["summary"]["flash_groups"], 2)

    def test_missing_data_and_zero_flow_accounted_for(self):
        model = heater()
        feed, product = model["streams"]
        product["component_flows_kmol_h"] = None
        model["streams"].append(dict(feed, id="OFF", molar_kmol_h=0., mass_kg_h=0.,
                                     component_flows_kmol_h={}, enthalpy_kw=0., entropy_kw_k=0.))
        plan = build_reference_plan(model)
        self.assertEqual(plan["already_covered"], {"OFF": "zero"})
        self.assertEqual(plan["unplannable"][0]["id"], "PRODUCT")
        self.assertEqual(plan["summary"]["pending_streams"], 2)
        self.assertEqual(plan["summary"]["plannable_streams"], 1)

    def test_representative_and_group_order_are_independent_of_input_order(self):
        model = heater()
        before = build_reference_plan(model)
        model["streams"].reverse()
        after = build_reference_plan(model)
        self.assertEqual(before, after)

    @unittest.skipUnless(ARCHIVE.exists(), "Requires the supplied case")
    def test_native_case_coverage_and_mixed_solid_work_instructions(self):
        model = load_case(ARCHIVE)
        plan = build_reference_plan(model)
        summary = plan["summary"]
        self.assertEqual((summary["total_streams"], summary["pending_streams"], summary["plannable_streams"]),
                         (171, 162, 162))
        self.assertEqual(summary["saved_ambient_candidates"], 0)
        ids = [sid for row in plan["streams"] for sid in row["member_ids"]]
        self.assertEqual(len(set(ids)), len(ids))
        self.assertEqual(len(ids), 162)
        self.assertEqual(len(plan["already_covered"]), 9)
        self.assertEqual(summary["flash_groups"], 76)
        water = next(row for row in plan["streams"] if "WAT" in row["member_ids"])
        self.assertNotIn("WAT-C", water["member_ids"])
        self.assertEqual(water["reference_option_hints"]["producer_free_water"], "YES")
        self.assertEqual(water["member_sources"]["WAT"]["producer_block"], "CONCEPTH")
        solid = next(row for row in plan["streams"] if "PROD4" in row["member_ids"])
        subflows = solid["substream_component_flows_kmol_h"]
        self.assertEqual(set(subflows), {"MIXED", "CISOLID"})
        self.assertGreater(subflows["CISOLID"]["C"], 0)


if __name__ == "__main__":
    unittest.main()
