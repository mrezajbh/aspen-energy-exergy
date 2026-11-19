import csv
import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from zipfile import ZipFile

from aspen_energy_exergy import analyze, load_case, plan_reference_flashes, write_flash_plan
from aspen_energy_exergy.cli import main


ROOT = Path(__file__).resolve().parents[1]
HEATER = ROOT / "examples/heater.json"
CASE = ROOT / "outputs/energy_exergy_aspen_20261111/source/Biogas to Jetfuel - FT.apwz"


class FlashPacketTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.model = load_case(HEATER)
        for stream in self.model["streams"]:
            stream["substreams"] = {"MIXED": {"component_flows_kmol_h": stream["component_flows_kmol_h"],
                                              "molar_kmol_h": stream["molar_kmol_h"]}}
        self.source = self.dir / "case.json"
        self.source.write_text(json.dumps(self.model), encoding="utf-8")
        self.plan = plan_reference_flashes(self.model)
        self.packet = write_flash_plan(self.plan, self.dir / "packet.zip")
        with ZipFile(self.packet) as archive:
            self.assertIsNone(archive.testzip())
            self.assertEqual(json.loads(archive.read("flash_plan.json")), self.plan)
            csv_in = csv.DictReader(io.StringIO(archive.read("flash_returns.csv").decode()))
            self.fields, self.rows = csv_in.fieldnames, list(csv_in)
        self.assertEqual(len(self.rows), 1)

    def write_csv(self, name, rows):
        path = self.dir / name
        with path.open("w", encoding="utf-8", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=self.fields)
            writer.writeheader()
            writer.writerows(rows)
        return path

    def completed(self):
        return {**self.rows[0], "T_K": 298.15, "P_Pa": 101325., "h0_J_mol": 0, "s0_J_molK": 0,
                "source_location": "Analytic hydrogen heater fixture", "aspen_reference_id": "TEST-REF",
                "representative_verified": "yes", "peers_verified": "yes"}

    def test_fill_one_flash_and_replan(self):
        blank = self.write_csv("blank.csv", self.rows)
        self.assertEqual(plan_reference_flashes(self.model, references=blank)["summary"]["flash_groups"], 1)
        filled = self.write_csv("filled.csv", [self.completed()])
        dest = self.dir / "results.zip"
        with redirect_stdout(io.StringIO()):
            self.assertEqual(main([str(self.source), "--references", str(blank), "--references", str(filled),
                                   "--native-splits-only", "-o", str(dest)]), 0)
        with ZipFile(dest) as archive:
            report = json.loads(archive.read("analysis.json"))
        self.assertEqual(report["coverage"]["split_quality"], {"matched": 2})
        self.assertEqual(report["streams"]["PRODUCT"]["reference_stream_id"], "FEED")
        s = report["streams"]["PRODUCT"]
        self.assertAlmostEqual(s["physical_exergy_kw"], s["enthalpy_kw"] - 298.15 * s["entropy_kw_k"])
        self.assertEqual(plan_reference_flashes(self.model, references=filled)["summary"]["flash_groups"], 0)
        with self.assertRaisesRegex(ValueError, "Multiple reference returns"):
            analyze(self.model, references=[filled, filled], estimate_splits=False)

    def test_unverified_peers_stale_inputs_and_partial_measurements_rejected(self):
        for updates, text in [
            ({"peers_verified": ""}, "peers_verified"), ({"P_Pa": 100000}, "reference P_Pa"),
            ({"operating_state_sha256": "0" * 64}, "operating state"),
            ({"representative_verified": ""}, "representative_verified"),
            ({"s0_J_molK": ""}, "s0_J_molK"), ({"h0_J_mol": "nan"}, "finite"),
            ({"applies_to": '["unknown"]'}, "unknown|invalid|Unknown"),
        ]:
            row = {**self.completed(), **updates}
            with self.subTest(updates=updates), self.assertRaisesRegex(ValueError, text):
                analyze(self.model, references=self.write_csv("bad.csv", [row]), estimate_splits=False)
        only_rep = {**self.completed(), "applies_to": "[]", "peers_verified": ""}
        report = analyze(self.model, references=self.write_csv("rep.csv", [only_rep]), estimate_splits=False)
        self.assertEqual(report["coverage"]["splits_available"], 1)

    def test_idempotent_packet_and_preserve_edited_json(self):
        self.assertEqual(write_flash_plan(self.plan, self.packet), self.packet)
        with ZipFile(self.packet) as archive:
            self.assertIn(self.plan["case_sha256"], archive.read("instructions.txt").decode())
            group = next(csv.DictReader(io.StringIO(archive.read("flash_groups.csv").decode())))
            self.assertEqual(set(json.loads(group["member_sources"])), set(json.loads(group["member_ids"])))
        with self.assertRaisesRegex(ValueError, "reference flash plan"):
            write_flash_plan({"streams": []}, self.dir / "invalid.json")
        path = write_flash_plan(self.plan, self.dir / "plan.json")
        path.write_text("operator changes", encoding="utf-8")
        with self.assertRaises(FileExistsError):
            write_flash_plan(self.plan, path)
        self.assertEqual(path.read_text(encoding="utf-8"), "operator changes")

    @unittest.skipUnless(CASE.exists(), "Local Aspen case required")
    def test_prior_native_sheet_matches_upgraded_native_reader(self):
        old = ROOT / "outputs/aspen_library/v0.2.0/Aspen_EXERGYFL_return.csv"
        plan = plan_reference_flashes(CASE, native_exergy=old)
        self.assertEqual(plan["summary"]["pending_streams"], 162)
        self.assertEqual(plan["summary"]["flash_groups"], 76)
