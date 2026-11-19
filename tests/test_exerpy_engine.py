"""Scientific integration checks for the real ExerPy arithmetic path."""
from copy import deepcopy
import importlib.abc
import sys
import unittest

from aspen_energy_exergy.engines.exerpy_engine import run_exerpy


def fixture():
    streams=[dict(id=s,mass_kg_h=3600,temperature_k=t,pressure_pa=101325,
                  enthalpy_kw=h,entropy_kw_k=.01,component_flows_kmol_h={'X':1})
             for s,t,h in [('IN',300,10),('OUT',350,20)]]
    model=dict(streams=streams,blocks=[dict(id='B',type='HEATER',inlets=['IN'],outlets=['OUT'])])
    properties={s:dict(total_exergy_kw=b,physical_exergy_kw=b-30,chemical_exergy_kw=30)
                for s,b in [('IN',50),('OUT',55)]}
    transfers=[dict(block='B',id='W',kind='work',rate_kw=10,exergy_kw=10),
               dict(block='B',id='Q',kind='heat',rate_kw=-2,exergy_kw=-1,basis='Supplied boundary heat exergy')]
    # The deliberately wrong old total proves the engine is recomputing from
    # primitive stream/transfer inputs rather than echoing the previous result.
    boundary=dict(id='B',inlets=['IN'],outlets=['OUT'],exergy_destruction_kw=777)
    balances=dict(blocks=[boundary],groups={},transfers=transfers,
                  plant={**boundary,'id':'plant','crossing_transfers':transfers})
    return model,properties,balances


class ExerPyIntegration(unittest.TestCase):
    def test_actual_ledger_signed_heat_and_perturbed_input(self):
        model,properties,balances=fixture()
        result=run_exerpy(model,properties,balances)
        self.assertEqual(result['boundaries']['B']['exergy_destruction_kw'],4)
        self.assertEqual(result['boundaries']['B']['difference_kw'],-773)
        modified=deepcopy(properties)
        modified['IN']['physical_exergy_kw']+=1
        modified['IN']['total_exergy_kw']+=1
        new=run_exerpy(model,modified,balances)
        self.assertEqual(new['boundaries']['B']['exergy_destruction_kw'],5)

    def test_unknown_split_keeps_total_without_fake_physical(self):
        model,properties,balances=fixture()
        properties['IN']['physical_exergy_kw']=None
        properties['IN']['chemical_exergy_kw']=None
        result=run_exerpy(model,properties,balances)
        row=result['stream_results']['IN']
        self.assertEqual(row['status'],'imported_total_no_split')
        self.assertIsNone(row['physical_exergy_kw'])
        self.assertIsNone(row['chemical_exergy_kw'])
        self.assertEqual(result['boundaries']['B']['exergy_destruction_kw'],4)

    def test_heat_and_work_perturbations_and_negative_balance(self):
        for index in (0,1):
            model,properties,balances=fixture()
            balances['transfers'][index]['exergy_kw']+=1
            result=run_exerpy(model,properties,balances)
            self.assertEqual(result['boundaries']['B']['exergy_destruction_kw'],5)
        model,properties,balances=fixture()
        properties['OUT']['total_exergy_kw']+=10
        properties['OUT']['physical_exergy_kw']+=10
        result=run_exerpy(model,properties,balances)
        self.assertEqual(result['boundaries']['B']['exergy_destruction_kw'],-6)
        self.assertTrue(any('negative' in w for w in result['boundaries']['B']['warnings']))

    def test_unknown_heat_is_not_native_validation(self):
        model,properties,balances=fixture()
        model['blocks'][0]['type']='COMPR'
        balances['transfers'][1]['exergy_kw']=None
        result=run_exerpy(model,properties,balances)
        self.assertEqual(result['boundaries']['B']['status'],'missing_input')
        self.assertFalse(result['component_checks']['B']['comparable_to_full_boundary'])
        self.assertTrue(any('unknown' in s for s in result['component_checks']['B']['limitations']))

    def test_no_coolprop_import_attempt(self):
        class Forbid(importlib.abc.MetaPathFinder):
            def find_spec(self,fullname,path=None,target=None):
                if fullname.lower().startswith('coolprop'):
                    raise AssertionError('CoolProp import attempted')
        guard=Forbid();sys.meta_path.insert(0,guard)
        try:
            result=run_exerpy(*fixture())
            self.assertTrue(result['provenance']['no_coolprop_loaded'])
            self.assertFalse(any(n.lower().startswith('coolprop') for n in sys.modules))
        finally:sys.meta_path.remove(guard)


if __name__=='__main__':unittest.main()
