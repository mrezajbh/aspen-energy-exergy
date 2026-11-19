import unittest
from unittest.mock import patch

from aspen_energy_exergy.performance import fuel_performance, _standard_lhv, _disable_optional_coolprop
from aspen_energy_exergy.engines import no_coolprop


def stream(name, species, enthalpy=0):
    return dict(id=name, component_flows_kmol_h=species, molar_kmol_h=sum(species.values()), enthalpy_kw=enthalpy)


class FuelPerformanceTest(unittest.TestCase):
    def test_reference_lookup_disables_optional_backend_without_import_probe(self):
        from thermo import coolprop
        _standard_lhv.cache_clear()
        # None is upstream's untested/initial state. Force it so a previous
        # test's availability cache cannot hide an unwanted import attempt.
        with patch.object(coolprop, '_has_CoolProp', None), no_coolprop() as attempts:
            methane, method = _standard_lhv('74-82-8')
            heavy, heavy_method = _standard_lhv('629-92-5')
            self.assertFalse(coolprop.has_CoolProp())
        self.assertEqual(attempts, [])
        self.assertEqual(methane, 802567.008)
        self.assertAlmostEqual(heavy, 11.7e6, delta=1e5)
        self.assertIn('ChEDL', method)
        self.assertIn('Estimate', heavy_method)

    def test_unknown_reference_backend_contract_fails_closed(self):
        from thermo import coolprop
        with patch.object(coolprop, 'has_CoolProp', None):
            with self.assertRaisesRegex(RuntimeError, 'availability contract'):
                _disable_optional_coolprop()

    def setUp(self):
        self.model = dict(
            components=[
                dict(id="water", cas="7732-18-5"), dict(id="fuel", cas="74-82-8"),
                dict(id="solid", name="carbon-graphite", atoms={"C": 1}),
                dict(id="unknown", name="unresolved pseudo oil"),
            ],
            streams=[
                stream("gas", {"water": 1, "fuel": 1}), stream("out", {"fuel": .9, "water": .2}),
                stream("steam", {"water": 3}, 0), stream("condensate", {"water": 3}, -10),
                stream("graphite", {"solid": 1}),
            ],
        )
        self.plant = dict(inlets=["gas", "steam", "graphite"], outlets=["out", "condensate"],
                          product_streams=["out"], heat_input_kw=2.0, work_input_kw=3.0)

    def test_condensing_utility_and_wet_fuel(self):
        f = fuel_performance(self.model, self.plant, ["gas", "graphite"], [("steam", "condensate")])
        self.assertAlmostEqual(f["components"]["fuel"]["lhv_j_mol"], 802567, delta=20)
        self.assertAlmostEqual(f["components"]["solid"]["lhv_j_mol"], 393474, delta=20)
        self.assertEqual(f["components"]["water"]["lhv_j_mol"], 0.0)
        feed, product = f["feed_lhv_kw"], f["product_lhv_kw"]
        self.assertAlmostEqual(f["efficiency_fraction"], product / (feed + 15))
        self.assertEqual(f["material_utility_heat_input_kw"], 10)
        self.assertAlmostEqual(f["streams"]["out"]["lhv_kw"], .9 * f["streams"]["gas"]["lhv_kw"])

    def test_receiving_utility_is_not_a_credit(self):
        self.model["streams"][3]["enthalpy_kw"] = 7
        f = fuel_performance(self.model, self.plant, ["gas"], [("steam", "condensate")])
        self.assertEqual(f["utility_pairs"][0]["enthalpy_drop_kw"], -7)
        self.assertEqual(f["material_utility_heat_input_kw"], 0)
        self.assertAlmostEqual(f["energy_input_kw"], f["feed_lhv_kw"] + 5)

    def test_unresolved_fuel_or_unmatched_utility_propagates(self):
        self.model["streams"][0] = stream("gas", {"fuel": .9, "unknown": .1})
        self.model["streams"][3] = stream("condensate", {"water": 1})
        f = fuel_performance(self.model, self.plant, ["gas"], [("steam", "condensate")])
        self.assertIsNone(f["feed_lhv_kw"])
        self.assertIsNone(f["efficiency_fraction"])
        self.assertIsNone(f["material_utility_heat_input_kw"])
        self.assertIsNotNone(f["product_lhv_kw"])
        self.assertTrue(f["diagnostics"])

    def test_selector_reuse_rejected(self):
        with self.assertRaisesRegex(ValueError, "overlaps"):
            fuel_performance(self.model, self.plant, ["gas", "steam"], [("steam", "condensate")])
        with self.assertRaisesRegex(ValueError, "wrong boundary"):
            fuel_performance(self.model, self.plant, ["out"])

    def test_case_missing_lower_heat_uses_phase_correct_formation_estimate(self):
        self.model["components"].append(dict(id="C19", cas="629-92-5"))
        self.model["streams"][1] = stream("out", {"C19": 1})
        f = fuel_performance(self.model, self.plant, ["gas"])
        self.assertIn("Estimate", f["components"]["C19"]["method"])
        self.assertAlmostEqual(f["components"]["C19"]["lhv_j_mol"], 11.7e6, delta=1e5)


if __name__ == "__main__":
    unittest.main()
