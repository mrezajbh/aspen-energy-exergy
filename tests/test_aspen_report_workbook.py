"""Regression checks for the Aspen workbook containing native physical exergy."""

from math import fsum
from pathlib import Path
import unittest

from aspen_energy_exergy import analyze


_BUNDLED = Path(__file__).resolve().parents[1] / "outputs/aspen_consistency_20261119/source_files"
_SOURCE = _BUNDLED if _BUNDLED.is_dir() else Path.home() / "Desktop"
_CASE = _SOURCE / "Biogas to Jetfuel - FT.apwz"
_WORKBOOK = _SOURCE / "EXPORT FULL REPORT WORKBOOK WITH EXERGY FLOW - MJ.xlsx"
_REPORT = _SOURCE / "REPORT TTXT.txt"


@unittest.skipUnless(all(p.is_file() for p in (_CASE, _WORKBOOK, _REPORT)),
                     "The three November Aspen exports are required")
class ExergyReportWorkbook(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = analyze(_CASE, native_exergy=_WORKBOOK, report_text=_REPORT, estimate_splits=False)

    def test_every_stream_signed_and_full_solid_coverage(self):
        r = self.report
        self.assertEqual(r["coverage"]["split_quality"], {"native": 164, "zero": 7})
        self.assertEqual(r["coverage"]["block_exergy_balances_available"], 99)
        self.assertEqual(len(r["export_summary"]["negative_whole_kw"]), 4)
        self.assertAlmostEqual(r["streams"]["BIOGAS"]["physical_exergy_kw"], -.03389367082886)
        for sid in ("PROD4", "PROD5", "REAL-OUT", "C-FEED", "FINES"):
            s = r["streams"][sid]
            parts = s["native_exergy"]["components_kw"]
            self.assertIn("CISOLID", parts)
            self.assertAlmostEqual(fsum(parts.values()), s["physical_exergy_kw"], places=9)
        jet = r["streams"]["JETFUEL"]
        self.assertAlmostEqual(jet["physical_exergy_kw"], 37.57389916029)
        self.assertAlmostEqual(jet["chemical_exergy_kw"], 13550.56136999, places=5)
        self.assertEqual(r["streams"]["KO-W"]["split_quality"], "zero")
        self.assertIsNone(r["streams"]["KO-W"]["export_evidence"]["native_whole_cal_s"])

    def test_report_stays_on_workbook_state(self):
        q = self.report["export_summary"]
        self.assertTrue(q["topology_checked"] and q["text_report_checked"])
        self.assertEqual(q["mixed_reported"], 162)
        self.assertEqual(q["solid_reported"], 5)
        self.assertLess(q["whole_minus_substreams_max_abs_kw"], 5e-12)
        self.assertEqual(self.report["source"]["material_state_source"], _WORKBOOK.name)
        self.assertAlmostEqual(self.report["plant"]["energy_residual_kw"], .02333749, places=6)
