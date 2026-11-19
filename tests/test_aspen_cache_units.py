"""Prevent silent SI interpretation of native COM display-unit caches."""

from pathlib import Path
import unittest

from aspen_energy_exergy.api import load_case
from aspen_energy_exergy.aspen import assert_cache_units


_ROOT = Path(__file__).resolve().parents[1]
_SOURCE = _ROOT / "outputs" / "aspen_consistency_20261119" / "source_files"
_OLD = _SOURCE / "Biogas to Jetfuel - FT.apwz"
_FRESH = _SOURCE / "native_rerun" / "native_case.bkp"
_NORMALIZED = _SOURCE / "native_rerun" / "model.json"


class CacheUnitDetection(unittest.TestCase):
    def test_si_column_is_accepted_and_celsius_cache_is_rejected(self):
        model = {
            "streams": [{"id": "top", "temperature_k": 278.15},
                        {"id": "bottom", "temperature_k": 502.715052}],
            "blocks": [{"id": "column", "type": "RADFRAC", "outlets": ["top", "bottom"],
                        "transfers": [{"id": "condenser", "process_temperature_k": 278.15},
                                      {"id": "reboiler", "process_temperature_k": 502.715052}]}],
        }
        assert_cache_units(model)
        for transfer in model["blocks"][0]["transfers"]:
            transfer["process_temperature_k"] -= 273.15
        with self.assertRaisesRegex(ValueError, "explicit native COM/report units"):
            assert_cache_units(model)

    @unittest.skipUnless(_OLD.is_file(), "Original native case fixture unavailable")
    def test_original_si_case_is_accepted(self):
        self.assertEqual(len(load_case(_OLD)["streams"]), 171)

    @unittest.skipUnless(_FRESH.is_file(), "Fresh native cache fixture unavailable")
    def test_fresh_display_unit_case_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "inconsistent display units"):
            load_case(_FRESH)

    @unittest.skipUnless(_NORMALIZED.is_file(), "Normalized native case fixture unavailable")
    def test_normalized_case_retains_explicit_native_units(self):
        blocks = {b["id"]: b for b in load_case(_NORMALIZED)["blocks"]}
        column = {t["id"]: t for t in blocks["T-601"]["transfers"]}
        self.assertAlmostEqual(column["condenser"]["rate_kw"], -127.036626)
        self.assertAlmostEqual(column["reboiler"]["rate_kw"], 42.5467275)
        self.assertAlmostEqual(column["condenser"]["process_temperature_k"], 278.15)
        self.assertAlmostEqual(column["reboiler"]["process_temperature_k"], 502.715052)
        self.assertAlmostEqual(blocks["HX-04"]["internal_heat_transfer_kw"], 632.599266)

