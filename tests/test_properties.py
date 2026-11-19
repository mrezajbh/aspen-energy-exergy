"""Datum, throughput and physical-property regression tests."""
import math
import unittest
from copy import deepcopy
from unittest.mock import patch

from aspen_energy_exergy._estimate import _fluid_delta, _molar_availability
from aspen_energy_exergy._exchange import case_hashes, digest
from aspen_energy_exergy.properties import REFERENCE_PLAN_SCHEMA, calculate_streams, equivalent_reference_streams

DIGEST = "abc1" * 16


def sample(name="METHANE", atoms=None):
    return {
        "metadata": {"source_sha256": DIGEST, "enthalpy_entropy_basis": "aspen_formation"},
        "components": [{"id": "X", "name": name, "atoms": atoms or {"C": 1, "H": 4}}],
        "streams": [{"id": "S", "property_method": "RKS-BM", "molar_kmol_h": 3.6,
                     "enthalpy_kw": -74.8, "entropy_kw_k": -0.08,
                     "component_flows_kmol_h": {"X": 3.6},
                     "temperature_k": 403.15, "pressure_pa": 2e5}],
    }


def flash_doc():
    return {"case_sha256": DIGEST,
            "streams": [{"id": "S", "component_ids": ["X"], "component_flows_kmol_h": [7.2],
                         "aspen_reference_id": "S-REF", "source_location": "ambient flash report",
                         "T_K": 298.15, "P_Pa": 101325., "property_method": "RKS-BM",
                         "quantity_basis": "flow_SI", "H0_kW": -150., "S0_kW_K": -0.17,
                         "same_aspen_property_options": True, "phase_convention": "equilibrium"}]}


def native_rate(value, scope="whole_stream"):
    return dict(physical_exergy_kw=value, scope=scope, export_source="Aspen stream report, row 6",
                case_sha256=DIGEST, reference_temperature_k=298.15, reference_pressure_pa=101325.)


def group_case():
    m = sample()
    a = m["streams"][0]
    a.update(mass_kg_h=57.6, substreams={"MIXED": {"component_flows_kmol_h": {"X": 3.6}}})
    b = dict(a, id="B", molar_kmol_h=7.2, mass_kg_h=115.2, component_flows_kmol_h={"X": 7.2},
             enthalpy_kw=-140., entropy_kw_k=-.12,
             substreams={"MIXED": {"component_flows_kmol_h": {"X": 7.2}}})
    m["streams"].append(b)
    ref = flash_doc()
    ref["streams"][0].update(applies_to=["B"], same_aspen_property_options_for_all=True)
    return m, ref


class StreamCalculations(unittest.TestCase):
    def test_shared_flash_scales_all_quantity_bases_and_preserves_provenance(self):
        for basis, values in (
                ("flow_SI", {"H0_kW": -150., "S0_kW_K": -.17}),
                ("molar_SI", {"h0_J_mol": -75000., "s0_J_molK": -85.}),
                ("specific_SI", {"h0_J_kg": -4687500., "s0_J_kgK": -5312.5})):
            with self.subTest(basis=basis):
                m, ref = group_case()
                r = ref["streams"][0]
                del r["H0_kW"], r["S0_kW_K"]
                r.update(quantity_basis=basis, **values)
                rows = calculate_streams(m, references=ref, estimate_splits=False)
                peer = rows["B"]
                self.assertEqual(peer["split_method"], "aspen_shared_ambient_flash")
                self.assertEqual(peer["split_quality"], "matched")
                self.assertEqual(peer["reference_stream_id"], "S")
                self.assertEqual(peer["aspen_reference_id"], "S-REF")
                self.assertEqual(peer["reference_source_location"], "ambient flash report")
                self.assertAlmostEqual(peer["physical_exergy_kw"], 10. - 298.15 * .05)
                self.assertAlmostEqual(peer["chemical_exergy_kw"], 2 * rows["S"]["chemical_exergy_kw"])
                native = {"B": native_rate(peer["physical_exergy_kw"])}
                direct = calculate_streams(m, references=ref, native_exergy=native, estimate_splits=False)["B"]
                self.assertEqual(direct["split_method"], "aspen_exergyfl")
                self.assertEqual(direct["reference_stream_id"], "S")

    def test_shared_flash_rejects_ambiguous_and_incompatible_peers(self):
        m, ref = group_case()
        a, b = m["streams"]
        record = ref["streams"][0]
        for ids, flag in ((["B"], False), (["B"], "true"), (["S"], True),
                          (["B","B"], True), ("B", True), (["missing"], True)):
            with self.subTest(ids=ids, flag=flag), self.assertRaises(ValueError):
                calculate_streams(m, references={**ref, "streams": [
                    {**record, "applies_to": ids, "same_aspen_property_options_for_all": flag}]})
        b["reference_options_key"] = "explicit free-water flash"
        with self.assertRaisesRegex(ValueError, "property options"):
            calculate_streams(m, references=ref)
        del b["reference_options_key"]
        saved = b.pop("substreams")
        with self.assertRaisesRegex(ValueError, "substreams"):
            calculate_streams(m, references=ref)
        b["substreams"] = saved
        b["mass_kg_h"] *= 1.01
        del record["H0_kW"], record["S0_kW_K"]
        record.update(quantity_basis="specific_SI", h0_J_kg=0., s0_J_kgK=0.)
        with self.assertRaisesRegex(ValueError, "molar mass"):
            calculate_streams(m, references=ref)
        self.assertTrue(equivalent_reference_streams(a,b))

    def test_same_bulk_composition_does_not_imply_same_substream_partition(self):
        m, ref = group_case()
        a,b = m["streams"]
        a["substreams"] = {"MIXED": {"component_flows_kmol_h": {"X": 3.5}},
                           "CISOLID": {"component_flows_kmol_h": {"X": .1}}}
        b["substreams"] = {"MIXED": {"component_flows_kmol_h": {"X": 6.}},
                           "CISOLID": {"component_flows_kmol_h": {"X": 1.2}}}
        self.assertFalse(equivalent_reference_streams(a,b))
        with self.assertRaisesRegex(ValueError, "substreams"):
            calculate_streams(m, references=ref)
        rows = calculate_streams(m)
        self.assertFalse(any("reference_anchor_id" in row for row in rows.values()))

    def test_template_ignores_blank_records_and_validates_incomplete_returns(self):
        m, ref = group_case()
        ref["plan_schema"] = REFERENCE_PLAN_SCHEMA
        blank = ref["streams"][0]
        blank.update(H0_kW=None, S0_kW_K=None, aspen_reference_id=None, source_location=None,
                     same_aspen_property_options_for_all=False, same_aspen_property_options=False)
        rows = calculate_streams(m, references=ref, estimate_splits=False)
        self.assertEqual([r["physical_exergy_kw"] for r in rows.values()], [None,None])
        partial = deepcopy(ref)
        partial["streams"][0]["H0_kW"] = 0.
        with self.assertRaises(ValueError):
            calculate_streams(m, references=partial)
        with self.assertRaisesRegex(ValueError, "plan schema"):
            calculate_streams(m, references={**ref, "plan_schema": "unknown"})

    def test_filled_plan_verifies_normalized_identity_and_substream_allocation(self):
        m, ref = group_case()
        del m["metadata"]["source_sha256"]
        ref.update(case_sha256=case_hashes(m)[0], plan_schema=REFERENCE_PLAN_SCHEMA)
        record = ref["streams"][0]
        record["operating_state_sha256"] = digest(m["streams"][0])
        record["substream_component_flows_kmol_h"] = {"MIXED": {"X": 1.}}
        rows = calculate_streams(m, references=ref, estimate_splits=False)
        self.assertEqual(rows["B"]["split_method"], "aspen_shared_ambient_flash")
        record["substream_component_flows_kmol_h"]["MIXED"] = {"Z": 1.}
        with self.assertRaisesRegex(ValueError, "substream component allocation"):
            calculate_streams(m, references=ref)
        record["operating_state_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "fingerprint"):
            calculate_streams(m, references=ref)
    def test_free_water_option_distinguishes_estimated_anchors(self):
        m, _ = group_case()
        m["streams"][0]["reference_options_key"] = "explicit free-water"
        rows = calculate_streams(m)
        self.assertFalse(any("reference_anchor_id" in row for row in rows.values()))

    def test_native_only_keeps_reported_physical_when_total_is_unavailable(self):
        m = sample()
        m["metadata"]["enthalpy_entropy_basis"] = "unverified"
        m["streams"][0].update(enthalpy_kw=None, component_flows_kmol_h=None)
        rec = native_rate(-0.01)
        row = calculate_streams(m, native_exergy={"S": rec}, estimate_splits=False)["S"]
        self.assertEqual(row["split_quality"], "native")
        self.assertEqual(row["physical_exergy_kw"], -0.01)
        self.assertIsNone(row["chemical_exergy_kw"])
        self.assertIsNone(row["total_exergy_kw"])
        self.assertEqual(row["native_exergy"], rec)

    def test_native_and_flash_must_agree(self):
        m = sample()
        physical = calculate_streams(m, references=flash_doc())["S"]["physical_exergy_kw"]
        row = calculate_streams(m, references=flash_doc(), native_exergy={"S": native_rate(physical + 1e-5)})["S"]
        self.assertEqual(row["split_method"], "aspen_exergyfl")
        self.assertEqual(row["physical_exergy_kw"], physical + 1e-5)
        self.assertAlmostEqual(row["chemical_exergy_kw"] + row["physical_exergy_kw"], row["total_exergy_kw"])
        with self.assertRaisesRegex(ValueError, "S: Aspen EXERGYFL.*disagrees"):
            calculate_streams(m, references=flash_doc(), native_exergy={"S": native_rate(physical + 2.)})
        with self.assertRaisesRegex(ValueError, "known stream IDs"):
            calculate_streams(m, native_exergy={"missing": native_rate(1.)})

    def test_native_graphite_scope_has_distinct_provenance(self):
        rec = native_rate(2., scope="mixed_with_saved_graphite")
        m = sample()
        m["components"].append({"id": "C", "name": "graphite", "atoms": {"C": 1}})
        m["streams"][0]["component_flows_kmol_h"] = {"X": 3.5, "C": .1}
        row = calculate_streams(m, native_exergy={"S": rec})["S"]
        self.assertEqual((row["split_method"], row["split_quality"]), ("aspen_exergyfl_with_saved_graphite", "hybrid"))
        self.assertAlmostEqual(row["chemical_exergy_kw"], row["total_exergy_kw"] - 2.)

    def test_native_anchor_is_preserved_and_conflicting_anchors_are_visible(self):
        m = sample()
        s = m["streams"][0]
        m["streams"] += [dict(s, id="B", temperature_k=320), dict(s, id="C", temperature_k=310)]
        rows = calculate_streams(m, native_exergy={"B": native_rate(1.)})
        self.assertEqual(rows["B"]["physical_exergy_kw"], 1.)
        self.assertEqual(rows["C"]["reference_anchor_id"], "B")
        self.assertEqual(rows["C"]["split_quality"], "estimated")
        rows = calculate_streams(m, native_exergy={"B": native_rate(1.), "C": native_rate(2.)})
        self.assertEqual(rows["C"]["physical_exergy_kw"], 2.)
        self.assertEqual(rows["B"]["physical_exergy_kw"], 1.)
        self.assertEqual(rows["S"]["split_method"], "thermo_peng_robinson")
        self.assertTrue(any("different molar chemical" in note for note in rows["S"]["notes"]))
        strict = calculate_streams(m, native_exergy={"B": native_rate(1.), "C": native_rate(2.)},
                                   estimate_splits=False)
        self.assertIsNone(strict["S"]["physical_exergy_kw"])
        differences = [r["reference_difference_kw"] for r in strict.values() if "reference_difference_kw" in r]
        self.assertEqual([abs(d) for d in differences], [1.])

    def test_native_flow_reference_rescales_and_closes(self):
        s = calculate_streams(sample(), references=flash_doc())["S"]
        env = .001 * (410260. + 4 * 118050.)
        self.assertAlmostEqual(s["total_exergy_kw"], -74.8 + 298.15 * .08 + env)
        self.assertAlmostEqual(s["physical_exergy_kw"], .2 - 298.15 * .005)
        self.assertAlmostEqual(s["chemical_exergy_kw"], s["total_exergy_kw"] - s["physical_exergy_kw"])
        self.assertEqual(s["split_quality"], "matched")
        self.assertEqual(s["element_flows_kmol_h"], {"C": 3.6, "H": 14.4})

    def test_reference_validation(self):
        for field, replacement in (("T_K", 300.), ("component_ids", ["other"]), ("property_method", "NRTL"),
                                   ("phase_convention", "liquid"), ("same_aspen_property_options", False),
                                   ("aspen_reference_id", ""), ("H0_kW", math.nan)):
            doc = flash_doc()
            doc["streams"][0][field] = replacement
            with self.subTest(field=field), self.assertRaises(ValueError):
                calculate_streams(sample(), references=doc)
        for field, replacement in (("case_sha256", "abcd" * 16), ("streams", flash_doc()["streams"] * 2)):
            doc = flash_doc()
            doc[field] = replacement
            with self.subTest(field=field), self.assertRaises(ValueError):
                calculate_streams(sample(), references=doc)

    def test_rejects_composition_mismatch(self):
        m = sample()
        m["components"].append({"id": "Y", "name": "NITROGEN", "atoms": {"N": 2}})
        doc = flash_doc()
        doc["streams"][0].update(component_ids=["X","Y"], component_flows_kmol_h=[3.6,3.6])
        with self.assertRaisesRegex(ValueError, "compositions differ"):
            calculate_streams(m, references=doc)

    def test_datum_and_disabling_estimates(self):
        m = sample()
        s = calculate_streams(m, estimate_splits=False)["S"]
        self.assertIsNone(s["physical_exergy_kw"])
        self.assertIsNotNone(s["total_exergy_kw"])
        m["metadata"]["enthalpy_entropy_basis"] = "unverified"
        s = calculate_streams(m, references=flash_doc())["S"]
        self.assertIsNone(s["chemical_exergy_kw"])
        self.assertIsNone(s["total_exergy_kw"])
        self.assertIsNotNone(s["physical_exergy_kw"])

    def test_near_dead_state_water_and_wet_steam(self):
        from iapws import IAPWS97
        m = sample("WATER", {"H": 2, "O": 1})
        stream = m["streams"][0]
        mass = 3.6 * 18.015268
        stream.update(mass_kg_h=mass, temperature_k=298.15, pressure_pa=101325., vapor_fraction=0)
        self.assertEqual(calculate_streams(m)["S"]["physical_exergy_kw"], 0)
        steam = IAPWS97(P=.4, x=.86)
        water = IAPWS97(P=.101325, T=298.15)
        stream.update(temperature_k=steam.T, pressure_pa=.4e6, vapor_fraction=.86)
        s = calculate_streams(m)["S"]
        self.assertEqual(s["split_method"], "iapws_if97")
        self.assertAlmostEqual(s["physical_exergy_kw"], mass / 3600 *
                               (steam.h - water.h - 298.15 * (steam.s - water.s)))
        self.assertAlmostEqual(s["total_exergy_kw"], s["physical_exergy_kw"] + s["chemical_exergy_kw"])

    def test_zero_and_graphite(self):
        m = sample("CARBON-GRAPHITE", {"C": 1})
        s = calculate_streams(m)["S"]
        self.assertEqual(s["split_quality"], "reference")
        self.assertAlmostEqual(s["chemical_exergy_kw"], 410.26)
        m["streams"][0].update(molar_kmol_h=0., component_flows_kmol_h={}, enthalpy_kw=None, entropy_kw_k=None)
        s = calculate_streams(m)["S"]
        self.assertEqual(s["split_quality"], "zero")
        self.assertEqual(s["total_exergy_kw"], 0)

    def test_bad_identification_and_unsupported_element(self):
        m = sample("WATER")  # inconsistent formula
        s = calculate_streams(m)["S"]
        self.assertEqual(s["split_quality"], "unavailable")
        self.assertIn("disagree", " ".join(s["notes"]))
        m = sample("ARGON", {"Ar": 1})
        s = calculate_streams(m, estimate_splits=False)["S"]
        self.assertIsNone(s["total_exergy_kw"])
        self.assertEqual(s["element_flows_kmol_h"], {"Ar": 3.6})

    def test_peng_robinson_ideal_gas_pressure_limit(self):
        # For gas at T=T0 and low pressure, specific exergy approaches RT ln(P/P0).
        m = sample("NITROGEN", {"N": 2})
        m["streams"][0].update(temperature_k=298.15, pressure_pa=2 * 101325.)
        s = calculate_streams(m)["S"]
        self.assertEqual(s["split_method"], "thermo_peng_robinson")
        self.assertAlmostEqual(s["physical_exergy_kw"], .001 * 8.31446 * 298.15 * math.log(2), delta=.02)
        self.assertEqual(s["split_quality"], "estimated")

    def test_same_composition_has_constant_molar_chemical_exergy(self):
        m = sample("WATER", {"H": 2, "O": 1})
        cold = m["streams"][0]
        cold.update(temperature_k=298.15, pressure_pa=101325., enthalpy_kw=-285.)
        hot = dict(cold, id="hot", molar_kmol_h=7.2, component_flows_kmol_h={"X": 7.2},
                   temperature_k=400., enthalpy_kw=-500., entropy_kw_k=-.12)
        m["streams"].append(hot)
        rows = calculate_streams(m)
        self.assertEqual(rows["hot"]["reference_anchor_id"], "S")
        self.assertAlmostEqual(rows["hot"]["chemical_exergy_kw"], 2 * rows["S"]["chemical_exergy_kw"])
        self.assertIsNotNone(rows["hot"]["independent_difference_kw"])
        self.assertAlmostEqual(rows["hot"]["physical_exergy_kw"],
                               (-500 + 298.15 * .12) - 2 * (-285 + 298.15 * .08))

    def test_failed_mixture_flash_uses_labeled_wilson_estimate(self):
        m = sample()
        m["components"].append({"id": "N", "name": "NITROGEN", "atoms": {"N": 2}})
        m["streams"][0]["component_flows_kmol_h"] = {"X": 1.8, "N": 1.8}
        _fluid_delta.cache_clear()
        self.addCleanup(_fluid_delta.cache_clear)
        with patch("aspen_energy_exergy._estimate._molar_availability", side_effect=Exception("flash failed")):
            s = calculate_streams(m)["S"]
        self.assertEqual(s["split_quality"], "estimated")
        self.assertEqual(s["split_method"], "thermo_pr_wilson_k")
        self.assertTrue(math.isfinite(s["physical_exergy_kw"]))

    def test_flash_cache_scales_flow_and_uses_exact_operating_state(self):
        _fluid_delta.cache_clear()
        _molar_availability.cache_clear()
        m = sample()
        first = calculate_streams(m)["S"]["physical_exergy_kw"]
        s = m["streams"][0]
        s.update(molar_kmol_h=7.2, component_flows_kmol_h={"X": 7.2})
        second = calculate_streams(m)["S"]["physical_exergy_kw"]
        self.assertAlmostEqual(second, 2 * first)
        self.assertEqual(_fluid_delta.cache_info().hits, 1)
        s["pressure_pa"] = math.nextafter(s["pressure_pa"], math.inf)
        calculate_streams(m)
        self.assertEqual(_fluid_delta.cache_info().misses, 2)
        self.assertEqual(_molar_availability.cache_info().hits, 1)  # common ambient state


if __name__ == "__main__":
    unittest.main()

