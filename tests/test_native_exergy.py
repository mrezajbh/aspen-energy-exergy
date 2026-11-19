import copy
import csv
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from aspen_energy_exergy import load_case
from aspen_energy_exergy.native_exergy import read_native_exergy, write_native_exergy_template


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "outputs/energy_exergy_aspen_20261111/source/Biogas to Jetfuel - FT.apwz"
LEGACY = ROOT / "outputs/aspen_library/v0.2.0/Aspen_EXERGYFL_return.csv"


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as source:
        return list(csv.DictReader(source))


def save(path, records):
    with path.open("w", encoding="utf-8-sig", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)


def completed(row, **changes):
    return {**row, "scope": "whole_stream", "whole_stream_EXERGYFL_kW": "-0.125",
            "export_source": "Run 2, property set PS-5", "reference_T_K": "298.15",
            "reference_P_Pa": "101325", "same_operating_state_confirmed": "yes", **changes}


class NativeExchange(unittest.TestCase):
    def setUp(self):
        workspace = TemporaryDirectory()
        self.addCleanup(workspace.cleanup)
        self.path = Path(workspace.name) / "native.csv"
        self.model = load_case(ROOT / "examples/heater.json")
        write_native_exergy_template(self.model, self.path)
        self.template = rows(self.path)

    def test_normalized_blank_partial_and_signed_whole_stream(self):
        self.assertEqual(read_native_exergy(self.path, self.model), {})
        self.assertEqual(self.template[0]["required_substreams"], "whole_stream")
        save(self.path, [completed(self.template[0]), self.template[1]])
        result = read_native_exergy(self.path, self.model)
        self.assertEqual(list(result), ["FEED"])
        self.assertEqual(result["FEED"]["physical_exergy_kw"], -0.125)
        self.assertEqual(result["FEED"]["export_source"], "Run 2, property set PS-5")
        self.assertEqual(result["FEED"]["components_kw"], {"whole_stream": -.125})
        reversed_model = {**self.model, "streams": list(reversed(self.model["streams"]))}
        self.assertEqual(result, read_native_exergy(self.path, reversed_model))

    def test_blank_template_generation_preserves_filled_file(self):
        self.assertEqual(write_native_exergy_template(self.model, self.path), self.path)
        save(self.path, [completed(self.template[0], whole_stream_EXERGYFL_kW="0")])
        before = self.path.read_bytes()
        with self.assertRaises(FileExistsError):
            write_native_exergy_template(self.model, self.path)
        self.assertEqual(self.path.read_bytes(), before)
        self.assertEqual(read_native_exergy(self.path, self.model)["FEED"]["physical_exergy_kw"], 0)

    def test_rejected_stale_or_ambiguous_values(self):
        for key, value in [
            ("property_method", "NRTL"), ("original_case_sha256", "e" * 64),
            ("operating_state_sha256", ""), ("operating_T_C", "26.1"), ("from_block", "other"),
            ("reference_T_K", "293.15"), ("reference_P_Pa", "100000"),
            ("whole_stream_EXERGYFL_kW", "NaN"), ("export_source", ""), ("MIXED_EXERGYFL_kW", "0"),
            ("same_operating_state_confirmed", "no"), ("id", "unknown"),
        ]:
            with self.subTest(key=key):
                save(self.path, [completed(self.template[0], **{key: value})])
                with self.assertRaises(ValueError):
                    read_native_exergy(self.path, self.model)
        save(self.path, [self.template[0], completed(self.template[0])])
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            read_native_exergy(self.path, self.model)
        save(self.path, [completed(self.template[0], scope="substreams", whole_stream_EXERGYFL_kW="",
                                   MIXED_EXERGYFL_kW="1")])
        with self.assertRaisesRegex(ValueError, "substream inventory"):
            read_native_exergy(self.path, self.model)

    def test_explicitly_wrong_units_are_rejected(self):
        row = completed(self.template[0])
        row["whole_stream_EXERGYFL_MW"] = row.pop("whole_stream_EXERGYFL_kW")
        save(self.path, [row])
        with self.assertRaisesRegex(ValueError, "rates in kW"):
            read_native_exergy(self.path, self.model)

    def test_composition_edit_under_preserved_source_hash(self):
        model = copy.deepcopy(self.model)
        model["metadata"]["source_sha256"] = "b" * 64
        original = self.path.with_name("source.csv")
        write_native_exergy_template(model, original)
        model["streams"][0]["component_flows_kmol_h"]["H2"] = 2.0
        with self.assertRaisesRegex(ValueError, "fingerprint"):
            read_native_exergy(original, model)


@unittest.skipUnless(ARCHIVE.exists() and LEGACY.exists(), "Saved Aspen archive and original return CSV required")
class SavedCaseNativeExchange(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = load_case(ARCHIVE)

    def test_original_sheet_and_new_backup_hash(self):
        self.assertEqual(read_native_exergy(LEGACY, self.model), {})
        with TemporaryDirectory() as d:
            p = write_native_exergy_template(self.model, Path(d) / "export.csv")
            requested = rows(p)
            self.assertEqual({r["id"] for r in requested}, {r["id"] for r in rows(LEGACY)})
            self.assertEqual(len(requested), 162)
            self.assertEqual(requested[0]["original_case_sha256"], self.model["metadata"]["backup_sha256"])
            backup = ARCHIVE.parent.parent / "native_backup_readonly.bkp"
            if backup.is_file():
                self.assertEqual(read_native_exergy(p, load_case(backup)), {})

    def test_graphite_and_stream_coverage(self):
        r = {row["id"]: row for row in rows(LEGACY)}
        with TemporaryDirectory() as d:
            out = Path(d) / "out.csv"
            a = completed(r["PROD4"], scope="mixed_with_saved_graphite", whole_stream_EXERGYFL_kW="",
                          MIXED_EXERGYFL_kW="15.25")
            b = completed(r["PROD5"], scope="substreams", whole_stream_EXERGYFL_kW="",
                          MIXED_EXERGYFL_kW="1", CISOLID_EXERGYFL_kW="0.03")
            c = completed(r["JETFUEL"], scope="substreams", whole_stream_EXERGYFL_kW="",
                          MIXED_EXERGYFL_kW="50", CISOLID_EXERGYFL_kW="0")
            save(out, [a, b, c, r["WAT"]])
            filled = read_native_exergy(out, self.model)
            self.assertEqual(set(filled), {"PROD4", "PROD5", "JETFUEL"})
            self.assertAlmostEqual(filled["PROD4"]["physical_exergy_kw"], 15.25 + 6.86334608474632)
            self.assertAlmostEqual(filled["PROD5"]["physical_exergy_kw"], 1.03)
            self.assertEqual(filled["JETFUEL"]["physical_exergy_kw"], 50.0)
            for row in ({**a, "CISOLID_EXERGYFL_kW": "1"}, {**b, "CISOLID_EXERGYFL_kW": ""},
                        {**c, "CISOLID_EXERGYFL_kW": "1"}):
                with self.subTest(sid=row["id"]):
                    save(out, [row])
                    with self.assertRaises(ValueError):
                        read_native_exergy(out, self.model)

    def test_whole_return_after_datum_confirmation(self):
        unknown = {**self.model, "metadata": {**self.model["metadata"], "enthalpy_entropy_basis": "unverified"}}
        with TemporaryDirectory() as d:
            path = write_native_exergy_template(unknown, Path(d) / "exergy.csv")
            prod4 = next(row for row in rows(path) if row["id"] == "PROD4")
            self.assertEqual(prod4["saved_graphite_physical_kW"], "")
            save(path, [completed(prod4)])
            self.assertEqual(read_native_exergy(path, self.model)["PROD4"]["physical_exergy_kw"], -.125)
            save(path, [completed(prod4, scope="mixed_with_saved_graphite", whole_stream_EXERGYFL_kW="",
                                  MIXED_EXERGYFL_kW="42")])
            with self.assertRaisesRegex(ValueError, "verified native graphite"):
                read_native_exergy(path, unknown)
            self.assertAlmostEqual(read_native_exergy(path, self.model)["PROD4"]["physical_exergy_kw"],
                                   42 + 6.86334608474632)
