"""Reference-flash option evidence extracted from the source deck."""

import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

from aspen_energy_exergy.aspen import _reference_option_hints, read_aspen


ARCHIVE = (Path(__file__).resolve().parents[1] / "outputs/energy_exergy_aspen_20261111/source"
           / "Biogas to Jetfuel - FT.apwz")


class OptionEvidence(unittest.TestCase):
    @unittest.skipUnless(ARCHIVE.is_file(), "Workspace Aspen regression archive unavailable")
    def test_native_request_packet_bytes_are_reproducible_across_hash_seeds(self):
        packets = []
        with TemporaryDirectory() as directory:
            for seed in ("1", "2"):
                packet = Path(directory) / f"flash-plan-{seed}.zip"
                subprocess.run(
                    [sys.executable, "-m", "aspen_energy_exergy", str(ARCHIVE),
                     "--write-flash-plan", str(packet)],
                    env={**os.environ, "PYTHONHASHSEED": seed},
                    check=True, capture_output=True,
                )
                packets.append(packet.read_bytes())
        self.assertTrue(packets[0] == packets[1], "Native flash request packet bytes differ across hash seeds")

    def test_quoted_delimiters_do_not_inject_block_options(self):
        deck = ('? FLOWSHEET GLOBAL ? \\ BLOCK BLKID = "H?1" \\ ? PROPERTIES MAIN ? '
                '\\ GPROPERTIES GBASEOPSET = IDEAL GFREE-WATER = STEAMNBS \\ ? BLOCK HEATER "H?1" ? '
                '\\ DESCRIPTION DESCRIPTION = "text ? BLOCK HEATER TEST ? \\ BLOCK-OPTION BLKOPFREWAT = YES" \\ ?')
        rows = [{"id": "OUT", "from_block": "H?1"}, {"id": "IN", "from_block": None}]
        options = _reference_option_hints(deck, rows)
        self.assertEqual(options["GLOBAL"]["GFREE-WATER"], "STEAMNBS")
        hints = rows[0]["reference_option_hints"]
        self.assertEqual(hints["section"], "GLOBAL")
        self.assertIsNone(hints["producer_free_water"])
        self.assertIsNone(rows[1]["reference_option_hints"]["section"])
        self.assertNotEqual(rows[0]["reference_options_key"], rows[1]["reference_options_key"])

    @unittest.skipUnless(ARCHIVE.is_file(), "Workspace Aspen regression archive unavailable")
    def test_saved_plant_property_scope_and_side_specific_options(self):
        case = read_aspen(ARCHIVE)
        streams = {r["id"]: r for r in case["streams"]}
        hints = {sid: row["reference_option_hints"] for sid, row in streams.items()}
        self.assertEqual(hints["WAT"]["producer_free_water"], "YES")
        self.assertEqual(hints["WAT"]["producer_block"], "CONCEPTH")
        self.assertIsNone(hints["WAT-C"]["producer_free_water"])
        self.assertEqual(hints["WAT-C"]["producer_block"], "COOLER")
        self.assertNotEqual(streams["WAT"]["reference_options_key"], streams["WAT-C"]["reference_options_key"])
        self.assertEqual(streams["WAT"]["reference_options_key"], streams["LIQ-WAR"]["reference_options_key"])
        self.assertEqual(hints["NAPHTHA"]["producer_free_water"], "NO")
        self.assertEqual(hints["NAPHTHA"]["producer_phase_options"]["NO-PHASE"], "3")
        self.assertEqual(hints["CWR-01"]["producer_property_options"]["OPSET"], "STEAMNBS")
        self.assertEqual(hints["CWR-01"]["exchanger_side"], "C")
        self.assertEqual(hints["TO-MEMB"]["producer_property_options"], {})
        self.assertEqual(hints["TO-MEMB"]["exchanger_side"], "H")
        self.assertEqual(hints["LPS-HEAD"]["section"], "STEAM")
        self.assertEqual(case["metadata"]["section_property_options"]["STEAM"]["SOLU-WATER"], "3")
