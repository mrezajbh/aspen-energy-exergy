import csv
import io
import json
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from zipfile import ZipFile

from aspen_energy_exergy import analyze, load_case, write_native_exergy_template, write_results
from aspen_energy_exergy.cli import main


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/heater.json"
PRIOR = ROOT / "outputs/energy_exergy_aspen_20261111"
ARCHIVE = PRIOR / "source/Biogas to Jetfuel - FT.apwz"


class PublicApiTest(unittest.TestCase):
    def test_cli_partial_native_return_round_trip(self):
        with TemporaryDirectory() as workspace, redirect_stdout(io.StringIO()):
            destination = Path(workspace)
            template = destination / "return.csv"
            self.assertEqual(main([str(EXAMPLE), "--write-exergy-template", str(template)]), 0)
            with template.open(encoding="utf-8-sig", newline="") as f:
                reader = csv.DictReader(f)
                fields, rows = reader.fieldnames, list(reader)
            feed = next(row for row in rows if row["id"] == "FEED")
            feed.update(export_source="Synthetic heater fixture at the reference state",
                        same_operating_state_confirmed="yes", reference_T_K="298.15",
                        reference_P_Pa="101325", scope="whole_stream", whole_stream_EXERGYFL_kW="0")
            with template.open("w", encoding="utf-8", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            report_path = destination / "complete.zip"
            self.assertEqual(main([str(EXAMPLE), "--native-exergy", str(template), "--native-splits-only",
                                   "-o", str(report_path)]), 0)
            with ZipFile(report_path) as z:
                report = json.loads(z.read("analysis.json"))
                self.assertIn("quality.csv", z.namelist())
                table = csv.DictReader(io.StringIO(z.read("streams.csv").decode()))
                csv_feed = next(r for r in table if r["id"] == "FEED")
                self.assertEqual(csv_feed["native_export_source"], feed["export_source"])
                self.assertEqual(json.loads(csv_feed["native_components_kw"]), {"whole_stream": 0})
            self.assertEqual(report["native_exergy_import"]["accepted_streams"], 1)
            self.assertEqual(report["streams"]["FEED"]["split_quality"], "native")
            self.assertEqual(report["streams"]["FEED"]["physical_exergy_kw"], 0)
            self.assertIsNone(report["streams"]["PRODUCT"]["physical_exergy_kw"])
            with self.assertRaises(FileExistsError):
                write_native_exergy_template(EXAMPLE, template)

    def test_unknown_temperature_keeps_first_law(self):
        model = load_case(EXAMPLE)
        transfer = model["blocks"][0]["transfers"][0]
        del transfer["temperature_k"]
        incomplete = analyze(model, estimate_splits=False)
        self.assertAlmostEqual(incomplete["plant"]["energy_residual_kw"], 0)
        self.assertIsNone(incomplete["plant"]["exergy_destruction_kw"])
        self.assertIn("HEATER/steam", " ".join(incomplete["plant"]["diagnostics"]))
        settings = {"transfers": [{"block": "HEATER", "id": "steam", "temperature_k": 420,
                                   "basis": "Steam measurement"}]}
        complete = analyze(model, settings=settings, estimate_splits=False)
        self.assertGreater(complete["plant"]["exergy_destruction_kw"], 0)
        self.assertNotIn("temperature_k", transfer)

    def test_profile_hash_blocks_unrelated_case(self):
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            analyze(EXAMPLE, settings="biogas-ft")

    def test_atomic_report(self):
        result = analyze(EXAMPLE, estimate_splits=False)
        with TemporaryDirectory() as workspace:
            path = write_results(result, Path(workspace) / "result.zip")
            with ZipFile(path) as archive:
                self.assertIsNone(archive.testzip())
                self.assertEqual(result, json.loads(archive.read("analysis.json")))
                self.assertIn("exergy_destruction_kw", archive.read("blocks.csv").decode())
            original = path.read_bytes()
            result["extra"] = float("nan")
            with self.assertRaises(ValueError):
                write_results(result, path)
            self.assertEqual(original, path.read_bytes())


@unittest.skipUnless(ARCHIVE.exists(), "Supplied Aspen archive and benchmark tables are required")
class SuppliedCaseTest(unittest.TestCase):
    def test_edit_one_heat_contact_preserves_rest_of_case_profile(self):
        model = load_case(ARCHIVE)
        old = analyze(model, settings="biogas-ft", estimate_splits=False)
        settings = {"base_profile": "biogas-ft", "transfers": [
            {"block": "HRU", "id": "duty", "temperature_k": 320., "basis": "Test heat-delivery scenario"}
        ]}
        new = analyze(model, settings=settings, estimate_splits=False)
        duty = next(t for t in old["transfers"] if t["block"] == "HRU" and t["id"] == "duty")
        delta = duty["rate_kw"] * 298.15 * (1 / duty["temperature_k"] - 1 / 320.)
        self.assertAlmostEqual(new["plant"]["exergy_destruction_kw"] - old["plant"]["exergy_destruction_kw"],
                               delta, places=8)
        for field in ("work_input_kw", "heat_input_kw", "energy_residual_kw", "product_streams"):
            self.assertEqual(old["plant"][field], new["plant"][field])
        self.assertEqual(new["groups"].keys(), old["groups"].keys())
        self.assertEqual(duty["temperature_k"], 308.15)
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            analyze(EXAMPLE, settings=settings, estimate_splits=False)
    def test_native_results_against_saved_analysis(self):
        result = analyze(ARCHIVE, estimate_splits=False)
        self.assertIn("matched by native backup SHA-256", result["settings_selection"])
        coverage = result["coverage"]
        self.assertEqual(coverage["streams"], 171)
        self.assertEqual(coverage["splits_available"], 9)
        self.assertEqual(coverage["block_exergy_balances_available"], 99)
        baseline = json.loads((PRIOR / "aspen_stream_results.json").read_text(encoding="utf-8"))
        for row in baseline:
            current = result["streams"][row["id"]]["total_exergy_kw"]
            self.assertAlmostEqual(current, row["total_exergy_kW"], delta=1e-7)
        plant = result["plant"]
        self.assertAlmostEqual(plant["exergy_destruction_kw"], 59574.00431962044, places=5)
        self.assertAlmostEqual(plant["energy_residual_kw"], 0.03595259635, places=6)
        self.assertAlmostEqual(plant["exergy_additivity_residual_kw"], 0., places=6)
        self.assertAlmostEqual(result["groups"]["cooler_water_section"]["exergy_destruction_kw"],
                               1693.5, delta=.1)
        self.assertEqual(sum(b["exergy_balance_status"] == "negative" for b in result["blocks"]), 8)
        self.assertAlmostEqual(result["fuel_performance"]["efficiency_fraction"], .1897504, delta=1e-4)
