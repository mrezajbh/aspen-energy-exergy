import csv
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from zipfile import ZipFile

from aspen_energy_exergy import analyze, load_case, write_results
from aspen_energy_exergy.quality import QUALITY_COLUMNS, build_quality
from aspen_energy_exergy.report import summary


HEATER = Path(__file__).resolve().parents[1] / "examples" / "heater.json"


def stream(sid, **fields):
    return {"id": sid, "total_exergy_kw": 1, "physical_exergy_kw": 1,
            "chemical_exergy_kw": 0, "split_quality": "estimated", **fields}


def empty_report():
    return {"streams": {}, "transfers": [], "blocks": [], "plant": {"id": "plant"},
            "source_warnings": [], "source": {"source_filename": "Example"}}


class QualityTests(unittest.TestCase):
    def test_same_composition_block_reports_signed_cross_option_difference(self):
        report = empty_report()
        for sid, free, chemical in [("W1", "YES", 20.), ("W2", None, 14.)]:
            report["streams"][sid] = stream(
                sid, chemical_exergy_kw=chemical, molar_kmol_h=10., component_flows_kmol_h={"W": 10.},
                property_method="RKS-BM", reference_options_key=sid,
                reference_option_hints={"producer_free_water": free},
                substreams={"MIXED": {"molar_kmol_h": 10., "component_flows_kmol_h": {"W": 10.}}})
        block = {"id": "HX", "inlets": ["W1"], "outlets": ["W2"]}
        report["blocks"] = [block]
        report["groups"] = {"also_HX": {**block, "id": "also_HX"}}
        findings = build_quality(report)
        self.assertEqual(len(findings), 1)
        row = findings[0]
        self.assertEqual((row["issue"], row["scope"], row["id"], row["value"]),
                         ("Cross-option chemical exergy difference", "block", "HX", -6.))
        self.assertIn("W1 to W2", row["evidence"])
        self.assertIn("YES to undeclared", row["evidence"])
        self.assertIn("option-matched", row["action"])
        last = report["streams"]["W2"]
        for fields in ({"reference_options_key": "W1"}, {"reference_options_key": None},
                       {"chemical_exergy_kw": 20.000001}, {"molar_kmol_h": 5.},
                       {"substreams": {"CISOLID": last["substreams"]["MIXED"]}}):
            with self.subTest(fields=fields):
                report["streams"]["W2"] = {**last, **fields}
                self.assertFalse(build_quality(report))

    def test_largest_model_differences_use_absolute_magnitude_and_preserve_sign(self):
        report = empty_report()
        for i, delta in enumerate([1, -6, 0, 4, -2, 5, 3]):
            sid = str(i)
            report["streams"][sid] = stream(sid, independent_difference_kw=delta, reference_anchor_id="ref")
        findings = build_quality(report)
        self.assertEqual([r["value"] for r in findings], [-6, 5, 4, 3, -2])
        self.assertTrue(all(r["unit"] == "kW" and "ref" in r["evidence"] for r in findings))

    def test_distinguishes_native_from_model_warnings(self):
        report = empty_report()
        for sid, quality in [("R", "estimated"), ("M", "native")]:
            report["streams"][sid] = stream(sid, split_quality=quality, split_method="shared_ambient_reference",
                                            independent_split_method="thermo_pr_wilson_k",
                                            active_phases=["Liquid1", "Liquid2"])
        report["source_warnings"] = ["R: * WARNING saved by Aspen"]
        findings = build_quality(report)
        self.assertEqual({r["issue"] for r in findings},
                         {"Wilson phase approximation", "Two-liquid phase mismatch", "Saved Aspen warning"})
        self.assertEqual({r["id"] for r in findings}, {"R"})

    def test_signed_findings_identify_boundary_and_smallest_enclosing_group(self):
        report = empty_report()
        report["blocks"] = [{"id": "A", "exergy_balance_status": "negative", "exergy_destruction_kw": -2,
                             "t0_entropy_generation_kw": -2}]
        report["groups"] = {"subsystem": {"id": "subsystem", "members": ["A", "B"], "exergy_destruction_kw": .4}}
        findings = build_quality(report)
        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0]["value"], -2)
        self.assertIn("+0.400 kW", findings[0]["evidence"])

    def test_legacy_report_reexport_and_csv_round_trip(self):
        model = load_case(HEATER)
        del model["blocks"][0]["transfers"][0]["temperature_k"]
        report = analyze(model, estimate_splits=False)
        report.pop("quality", None)
        self.assertIn("Missing heat boundary temperature", summary(report))
        with TemporaryDirectory() as tmp:
            output = write_results(report, Path(tmp) / "case.zip")
            with ZipFile(output) as z:
                reader = csv.DictReader(io.StringIO(z.read("quality.csv").decode()))
                findings = list(reader)
                self.assertEqual(tuple(reader.fieldnames), QUALITY_COLUMNS)
                row = next(r for r in findings if r["issue"] == "Missing heat boundary temperature")
                self.assertEqual((row["scope"], row["id"], row["value"], row["unit"]),
                                 ("transfer", "HEATER/steam", "", "K"))
                self.assertEqual(len(findings), len(json.loads(z.read("analysis.json"))["quality"]))
        self.assertNotIn("quality", report)

    def test_complete_clean_stream_has_no_findings(self):
        report = empty_report()
        report["streams"]["R"] = stream("R", split_quality="native", active_phases=["Liquid1", "Liquid2"])
        self.assertEqual(build_quality(report), [])
