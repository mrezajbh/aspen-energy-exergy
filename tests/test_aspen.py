"""Regression checks for units, signed transfers and the saved Aspen file interface."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from zipfile import ZipFile

from aspen_energy_exergy.aspen import read_aspen


_CASE = Path(__file__).resolve().parents[1] / "outputs" / "energy_exergy_aspen_20261111"
_ARCHIVE = _CASE / "source" / "Biogas to Jetfuel - FT.apwz"


class InputErrors(unittest.TestCase):
    def test_version_and_missing_results(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / "case.bkp"
            path.write_bytes(b'MM "40.0"')
            with self.assertRaisesRegex(ValueError, "version"):
                read_aspen(path)
            path.write_bytes(b'MM "41.0"')
            with self.assertRaisesRegex(ValueError, "results"):
                read_aspen(path)

    def test_apwz_members(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / "case.apwz"
            with ZipFile(path, "w") as archive:
                archive.writestr("case.bkp", b'x')
                archive.writestr("other.bkp", b'x')
            with self.assertRaisesRegex(ValueError, "one .bkp"):
                read_aspen(path)


@unittest.skipUnless(_ARCHIVE.is_file(), "Workspace Aspen regression archive unavailable")
class SavedAspenCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.case = read_aspen(_ARCHIVE)
        cls.streams = {s["id"]: s for s in cls.case["streams"]}
        cls.blocks = {b["id"]: b for b in cls.case["blocks"]}

    def test_backup_and_archive_return_same_values(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / "renamed.bkp"
            with ZipFile(_ARCHIVE) as archive:
                path.write_bytes(archive.read(self.case["metadata"]["source_member"]))
            from_bkp = read_aspen(path)
            for key in ("components", "streams", "blocks", "energy_links", "warnings"):
                self.assertEqual(self.case[key], from_bkp[key], key)
            self.assertEqual(self.case["metadata"]["backup_sha256"], from_bkp["metadata"]["backup_sha256"])
            self.assertEqual(from_bkp["metadata"]["enthalpy_entropy_basis"], "aspen_formation")

    def test_inventory_and_formation_basis(self):
        self.assertEqual((len(self.streams), len(self.blocks), len(self.case["components"])), (171, 99, 65))
        self.assertEqual(len(self.case["energy_links"]), 12)
        self.assertEqual(len(self.case["warnings"]), 2)
        self.assertEqual(self.case["metadata"]["enthalpy_entropy_basis"], "aspen_formation")
        stream = self.streams["1"]
        self.assertAlmostEqual(stream["temperature_k"], 333.061420399015)
        self.assertEqual(stream["pressure_pa"], 150000.)
        self.assertAlmostEqual(stream["enthalpy_kw"], -9670.682297983136)
        self.assertAlmostEqual(stream["entropy_kw_k"], -10.215768247764416)
        self.assertAlmostEqual(self.streams["C-FEED"]["component_flows_kmol_h"]["C"], 46.50497206572708)
        self.assertEqual(self.streams["BTTM1"]["entropy_kw_k"], 0.)
        self.assertIsNone(self.streams["BTTM1"]["temperature_k"])

    def test_phase_and_substream_identity(self):
        self.assertEqual(self.case["components"][0]["id"], "H2O")
        self.assertEqual(self.streams["WAT"]["active_phases"], ["Liquid2"])
        self.assertEqual(self.streams["WAT-C"]["active_phases"], ["Liquid1"])
        self.assertEqual(self.streams["NAPHTHA"]["active_phases"], ["Liquid1", "Liquid2"])
        self.assertEqual(self.streams["COL-FEED"]["active_phases"], ["Vapor", "Liquid1", "Liquid2"])
        self.assertTrue(all(s["active_phases"] is not None for s in self.streams.values()))
        solid = self.streams["PROD4"]["substreams"]["CISOLID"]
        self.assertEqual(set(solid["component_flows_kmol_h"]), {"C"})
        self.assertAlmostEqual(solid["enthalpy_kw"] - 298.15 * solid["entropy_kw_k"], 6.86334608474632)

    def test_recovered_heat_and_shaft_signs(self):
        donor, receiver = (self.blocks[b]["transfers"] for b in ("R-FT", "COIL-WOU"))
        self.assertEqual(len(donor), 1)
        self.assertEqual(len(receiver), 1)
        self.assertEqual(donor[0]["link_id"], receiver[0]["link_id"])
        self.assertAlmostEqual(donor[0]["rate_kw"], -4851.4869586246905)
        self.assertEqual(donor[0]["rate_kw"], -receiver[0]["rate_kw"])
        self.assertEqual(self.blocks["HX-01"]["transfers"], [])
        self.assertGreater(self.blocks["HX-01"]["internal_heat_transfer_kw"], 0.)
        for bid, expected in (("K-101A", 263.171520527503), ("ST-TURB", -3262.23265171954)):
            power, loss = self.blocks[bid]["transfers"]
            self.assertEqual(power["kind"], "work")
            self.assertAlmostEqual(power["rate_kw"], expected)
            self.assertLess(loss["rate_kw"], 0.)
        cond, reb = self.blocks["T-601"]["transfers"]
        self.assertAlmostEqual(cond["rate_kw"], -127.036431)
        self.assertAlmostEqual(reb["rate_kw"], 42.5465051)
        self.assertIsNone(cond["temperature_k"])
        self.assertAlmostEqual(cond["process_temperature_k"], 278.15)
