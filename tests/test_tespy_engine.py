"""Tests of native TESPy invocation, missing inputs and independent perturbations."""
import copy
import unittest
from unittest.mock import patch

from aspen_energy_exergy.engines.tespy_engine import run_tespy,_native_runtime


class TespyImportedStateTest(unittest.TestCase):
    def setUp(self):
        self.model={'streams':[
            {'id':'in','enthalpy_kw':-100.,'mass_kg_h':3600.},
            {'id':'out','enthalpy_kw':-90.,'mass_kg_h':3600.}],
            'blocks':[{'id':'HX','type':'HEATER','inlets':['in'],'outlets':['out'],
                'transfers':[{'id':'Q','kind':'heat','rate_kw':10.}]}]}
        self.boundaries={'blocks':[{'id':'HX','inlets':['in'],'outlets':['out'],
            'energy_residual_kw':987654.}], 'groups':[],
            'plant':{'id':'plant','inlets':['in'],'outlets':['out'],
                     'crossing_transfers':[{'id':'Q','kind':'heat','rate_kw':10.}]}}

    def test_actual_native_equation_called_and_old_answer_ignored(self):
        with _native_runtime() as (native,guard):
            cls=native['HeatBus']
            original=cls.energy_balance_func
            calls=[]
            def observed(obj):
                calls.append(obj.label)
                return original(obj)
            with patch.object(cls,'energy_balance_func',observed):
                r=run_tespy(self.model,{},self.boundaries)
        self.assertEqual(len(calls),2)
        self.assertEqual(r['boundaries']['HX']['energy_residual_kw'],0.)
        self.assertEqual(r['component_checks']['HX']['energy_residual_kw'],0.)
        self.assertEqual(r['provenance']['CoolProp_import_attempts'],[])
        self.assertFalse(r['provenance']['CoolProp_loaded'])

    def test_output_enthalpy_perturbation_reaches_native_engine(self):
        self.model['streams'][1]['enthalpy_kw']=-89.
        r=run_tespy(self.model,{},self.boundaries)
        self.assertEqual(r['boundaries']['HX']['energy_residual_kw'],-1.)
        self.assertEqual(r['component_checks']['HX']['energy_residual_kw'],-1.)

    def test_missing_enthalpy_stays_unavailable(self):
        self.model['streams'][1]['enthalpy_kw']=None
        r=run_tespy(self.model,{},self.boundaries)
        self.assertIsNone(r['boundaries']['HX']['energy_residual_kw'])
        self.assertEqual(r['boundaries']['HX']['status'],'unavailable')

    def test_missing_transfer_inventory_is_not_assumed_adiabatic(self):
        del self.model['blocks'][0]['transfers']
        del self.boundaries['plant']['crossing_transfers']
        r=run_tespy(self.model,{},self.boundaries)
        self.assertTrue(all(v['energy_residual_kw'] is None for v in r['boundaries'].values()))

    def test_missing_mass_does_not_remove_valid_energy_rate_ledger(self):
        self.model['streams'][0]['mass_kg_h']=None
        r=run_tespy(self.model,{},self.boundaries)
        self.assertEqual(r['boundaries']['HX']['energy_residual_kw'],0.)
        self.assertEqual(r['component_checks']['HX']['status'],'not_applicable')


if __name__=='__main__':unittest.main()
