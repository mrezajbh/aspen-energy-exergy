import math
import unittest

from aspen_energy_exergy.balances import calculate_balances


def stream(name, enthalpy, entropy, mass=1.0):
    return dict(id=name, enthalpy_kw=enthalpy, entropy_kw_k=entropy, mass_kg_h=mass)


def block(name, inlets, outlets, transfers):
    return dict(id=name, inlets=inlets, outlets=outlets, transfers=transfers)


def heat(name, rate, temperature=None, link=None):
    return dict(id=name, kind="heat", rate_kw=rate, temperature_k=temperature, link_id=link)


class BalancesTest(unittest.TestCase):
    def setUp(self):
        self.states = [stream("feed", 20, .05), stream("mid", 14, .04), stream("product", 20, .06)]
        self.availability = {s: {"total_exergy_kw": x} for s, x in [("feed", 15), ("mid", 12), ("product", 12)]}
        self.blocks = [
            block("donor", ["feed"], ["mid"], [heat("coil", -6, 500, "Q")]),
            block("recipient", ["mid"], ["product"], [heat("coil", 6, 400, "Q")]),
        ]

    def run_balance(self, **kwargs):
        model = dict(streams=self.states, blocks=self.blocks)
        return calculate_balances(model, self.availability, t0_k=300, **kwargs)

    def test_internal_link_with_two_temperatures(self):
        result = self.run_balance(product_streams=["product"])
        plant = result["plant"]
        self.assertEqual(plant["internal_material_streams"], ["mid"])
        self.assertEqual(plant["internal_energy_links"], ["Q"])
        self.assertEqual(plant["energy_residual_kw"], 0)
        self.assertEqual(plant["heat_input_kw"], 0)
        self.assertAlmostEqual(plant["internal_link_exergy_destruction_kw"], .9)
        self.assertAlmostEqual(plant["block_exergy_destruction_sum_kw"], 2.1)
        self.assertAlmostEqual(plant["exergy_additivity_residual_kw"], 0)
        self.assertEqual(plant["exergy_destruction_kw"], 3)
        self.assertAlmostEqual(plant["t0_entropy_generation_kw"], 3)
        self.assertEqual(plant["product_exergy_fraction_of_input"], .8)

    def test_internal_missing_state_does_not_poison_plant(self):
        self.availability["mid"]["total_exergy_kw"] = None
        self.blocks[0]["transfers"][0]["temperature_k"] = None
        self.blocks[1]["transfers"][0]["rate_kw"] = None
        result = self.run_balance()
        self.assertIsNone(result["blocks"][0]["exergy_destruction_kw"])
        self.assertIsNone(result["blocks"][1]["energy_residual_kw"])
        self.assertEqual(result["plant"]["exergy_destruction_kw"], 3)
        self.assertEqual(result["plant"]["energy_balance_status"], "closed")
        self.assertIsNone(result["plant"]["block_exergy_destruction_sum_kw"])

    def test_crossing_and_enclosed_links_on_group_boundaries(self):
        self.blocks[0]["transfers"].append(dict(id="drive", kind="work", rate_kw=-1, link_id="W"))
        self.blocks[1]["transfers"].append(dict(id="drive", kind="work", rate_kw=1, link_id="W"))
        groups = self.run_balance(groups={"left": ["donor"], "right": ["recipient"],
                                          "both": ["donor", "recipient"]})["groups"]
        left, right, both = groups["left"], groups["right"], groups["both"]
        self.assertEqual(left["inlets"], ["feed"])
        self.assertEqual(left["outlets"], ["mid"])
        self.assertEqual(left["heat_rejected_kw"], 6)
        self.assertEqual(left["work_output_kw"], 1)
        self.assertAlmostEqual(left["heat_exergy_net_kw"], -2.4)
        self.assertAlmostEqual(right["heat_exergy_net_kw"], 1.5)
        self.assertEqual(len(right["crossing_transfers"]), 2)
        self.assertEqual(both["internal_energy_links"], ["Q", "W"])
        self.assertEqual(both["internal_material_streams"], ["mid"])
        self.assertEqual(both["crossing_transfers"], [])
        self.assertAlmostEqual(both["exergy_additivity_residual_kw"], 0)
        self.assertAlmostEqual(both["internal_link_exergy_destruction_kw"], .9)

    def test_malformed_group_names_and_members(self):
        for groups in ({"bad": ["donor", "donor"]}, {"bad": []}, {"bad": ["unknown"]},
                       {"bad": "donor"}, {"": ["donor"]}, [["bad", "donor"]]):
            with self.subTest(groups=groups), self.assertRaises(ValueError):
                self.run_balance(groups=groups)

    def test_subambient_heat_and_work_signs(self):
        self.blocks[0]["transfers"].extend([heat("cold", -10, 250), dict(id="shaft", kind="work", rate_kw=-2)])
        result = self.run_balance()
        plant = result["plant"]
        self.assertEqual(plant["heat_rejected_kw"], 10)
        self.assertEqual(plant["work_output_kw"], 2)
        self.assertAlmostEqual(plant["heat_exergy_net_kw"], 2)
        self.assertAlmostEqual(plant["exergy_input_kw"], 17)
        self.assertAlmostEqual(plant["exergy_output_kw"], 14)
        self.assertEqual(plant["energy_residual_kw"], -12)
        self.assertEqual(plant["energy_balance_status"], "residual")

    def test_unknown_external_temperature_only_propagates_to_dependent_quantities(self):
        self.blocks[0]["transfers"].append(heat("hot", 5, None))
        plant = self.run_balance()["plant"]
        self.assertEqual(plant["energy_residual_kw"], 5)
        self.assertIsNone(plant["exergy_destruction_kw"])
        self.assertIsNone(plant["entropy_generation_kw_k"])
        self.assertIn("Missing boundary temperature for donor/hot", plant["diagnostics"])

    def test_atom_balance_and_exergy_identity(self):
        for row in self.availability.values():
            row.update(element_flows_kmol_h={"C": 2.0}, element_reference_kw=10.0)
        self.availability["product"].update(total_exergy_kw=11.0, element_reference_kw=9.0,
                                            element_flows_kmol_h={"C": 1.8})
        self.blocks[0]["internal_heat_transfer_kw"] = 94.0
        result = self.run_balance()
        plant = result["plant"]
        self.assertEqual(result["blocks"][0]["internal_heat_transfer_kw"], 94.0)
        self.assertAlmostEqual(plant["element_residual_kmol_h"]["C"], .2)
        self.assertEqual(plant["element_balance_status"], "residual")
        self.assertEqual(plant["element_reference_residual_kw"], 1)
        self.assertAlmostEqual(plant["exergy_identity_residual_kw"], 0)
        self.assertEqual(plant["heat_input_kw"], 0)

    def test_zero_heat_needs_no_temperature(self):
        self.blocks[0]["transfers"].append(heat("zero", 0))
        self.assertEqual(self.run_balance()["plant"]["exergy_destruction_kw"], 3)

    def test_temperature_range_crossing_ambient_retains_both_exergy_directions(self):
        self.blocks[0]["transfers"].append(dict(id="external", kind="heat", rate_kw=10,
                                                temperature_range_k=[350., 250.]))
        result = self.run_balance()
        imported = next(t for t in result["transfers"] if t["id"] == "external")
        cold = 5 - 30 * math.log(300 / 250)
        hot = 5 - 30 * math.log(350 / 300)
        self.assertAlmostEqual(imported["temperature_k"], 100 / math.log(350 / 250))
        self.assertAlmostEqual(imported["entropy_kw_k"], 0.1 * math.log(350 / 250))
        self.assertAlmostEqual(imported["exergy_kw"], cold + hot)
        self.assertAlmostEqual(imported["exergy_input_kw"], hot)
        self.assertAlmostEqual(imported["exergy_output_kw"], -cold)
        self.assertAlmostEqual(result["plant"]["exergy_input_kw"], 15 + hot)
        self.assertAlmostEqual(result["plant"]["exergy_output_kw"], 12 - cold)
        self.assertAlmostEqual(result["plant"]["exergy_additivity_residual_kw"], 0)

    def test_isothermal_limit_and_invalid_heat_ranges(self):
        transfer = dict(id="external", kind="heat", rate_kw=-1, temperature_range_k=[400., 400.])
        self.blocks[0]["transfers"].append(transfer)
        actual = next(t for t in self.run_balance()["transfers"] if t["id"] == "external")
        self.assertAlmostEqual(actual["exergy_kw"], -0.25)
        transfer["temperature_range_k"][1] = math.nextafter(400., math.inf)
        actual = next(t for t in self.run_balance()["transfers"] if t["id"] == "external")
        self.assertAlmostEqual(actual["temperature_k"], 400.)
        for limits in ([0, 300], [300], [True, 320], [300, None]):
            transfer["temperature_range_k"] = limits
            with self.subTest(limits=limits), self.assertRaises(ValueError):
                self.run_balance()
        transfer.update(temperature_k=301., temperature_range_k=[300., 305.])
        with self.assertRaisesRegex(ValueError, "no scalar"):
            self.run_balance()

    def test_invalid_link_and_ambiguous_topology(self):
        self.blocks[0]["transfers"][0]["rate_kw"] = -4
        with self.assertRaisesRegex(ValueError, "do not cancel"):
            self.run_balance()
        self.blocks[1]["inlets"].append("feed")
        with self.assertRaisesRegex(ValueError, "Multiple inlets"):
            self.run_balance()

    def test_product_must_cross_plant_boundary(self):
        with self.assertRaisesRegex(ValueError, "distinct external material outlet"):
            self.run_balance(product_streams=["mid"])

    def test_preserve_negative_destruction(self):
        self.availability["product"]["total_exergy_kw"] = 16
        plant = self.run_balance()["plant"]
        self.assertEqual(plant["exergy_destruction_kw"], -1)
        self.assertEqual(plant["exergy_balance_status"], "negative")


if __name__ == "__main__":
    unittest.main()
