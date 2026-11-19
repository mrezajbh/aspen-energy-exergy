"""Public imported-state mode must execute and publish the real upstream engines."""
from copy import deepcopy
from pathlib import Path
import unittest
from unittest.mock import patch

from aspen_energy_exergy import analyze, load_case
from aspen_energy_exergy.engines import no_coolprop, promote_external_results


EXAMPLE = Path(__file__).resolve().parents[1] / 'examples/heater.json'


class ExternalApiTest(unittest.TestCase):
    def test_invalid_options_are_rejected_before_any_inputs_or_properties(self):
        for options in ({'balance_engine': 'unknown'}, {'balance_engine': []},
                        {'balance_engine': 'exerpy-tespy'},
                        {'balance_engine': 'exerpy-tespy', 'estimate_splits': 0}):
            with self.subTest(options=options), patch('aspen_energy_exergy.api._analysis_inputs') as inputs, \
                    patch('aspen_energy_exergy.properties.calculate_streams') as properties:
                with self.assertRaises(ValueError):
                    analyze('this-source-must-not-be-read', **options)
                inputs.assert_not_called()
                properties.assert_not_called()

    def test_real_engines_publish_stream_block_group_and_plant_results(self):
        with no_coolprop() as attempts:
            report = analyze(EXAMPLE, estimate_splits=False, balance_engine='exerpy-tespy',
                             settings={'groups': {'heater_package': ['HEATER']},
                                       'product_streams': ['PRODUCT']})
        external = report['external_engines']
        self.assertEqual(attempts, [])
        self.assertEqual(external['coolprop'], {'import_attempts': [], 'loaded': False})
        for record in [*report['blocks'], *report['groups'].values(), report['plant']]:
            for engine, field in [('ExerPy', 'exergy_destruction_kw'), ('TESPy', 'energy_residual_kw')]:
                native = external[engine]['boundaries'][record['id']]
                self.assertEqual(record[field], native[field])
                self.assertEqual(record[field + '_engine'], engine)
                self.assertEqual(record[field + '_provenance']['method'], native['method'])
        self.assertGreater(report['plant']['exergy_destruction_kw'], 0)
        self.assertAlmostEqual(report['plant']['energy_residual_kw'], 0)
        for sid, stream in report['streams'].items():
            native = external['ExerPy']['stream_results'][sid]
            self.assertEqual(stream['total_exergy_kw'], native['total_exergy_kw'])
            self.assertEqual(stream['total_exergy_aggregation'], native['method'])
        # The fixture has no supplied product split: retaining a total is not
        # mislabelled as native ExerPy physical-property generation.
        self.assertEqual(report['streams']['PRODUCT']['total_exergy_kw_provenance']['engine'], 'imported')
        self.assertIsNone(report['streams']['PRODUCT']['physical_exergy_kw'])
        self.assertEqual(report['balance_engine'], 'exerpy-tespy')
        self.assertIn('ExerPy', report['software_versions'])
        self.assertIn('TESPy', report['software_versions'])
        # The shared helper rejects disagreements before mutating the report.
        broken = deepcopy(report)
        broken['external_engines']['comparisons'][0]['agreement_within_tolerance'] = False
        previous = deepcopy(broken)
        with self.assertRaisesRegex(ValueError, 'External balance disagreement'):
            promote_external_results(broken)
        self.assertEqual(broken, previous)

    def test_missing_boundary_data_stays_unavailable(self):
        model = load_case(EXAMPLE)
        del model['blocks'][0]['transfers'][0]['temperature_k']
        report = analyze(model, estimate_splits=False, balance_engine='exerpy-tespy')
        self.assertIsNone(report['plant']['exergy_destruction_kw'])
        self.assertEqual(report['plant']['exergy_balance_status'], 'incomplete')
        self.assertEqual(report['plant']['exergy_destruction_kw_provenance']['status'], 'missing_input')
        self.assertAlmostEqual(report['plant']['energy_residual_kw'], 0)


if __name__ == '__main__':
    unittest.main()
