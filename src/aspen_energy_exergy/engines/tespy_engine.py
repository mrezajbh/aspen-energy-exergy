"""Evaluate TESPy's native energy equations on imported Aspen state records.

TESPy 0.12 delegates exergy assessment to ExerPy. This adapter uses the native
TESPy energy-bus residual for every boundary and, where port topology permits,
native equipment energy residuals. It does not call a property engine, solve a
TESPy process network, or claim independent reactor/separator physics validation.
"""
from __future__ import annotations

from contextlib import contextmanager
import importlib.abc
import json
from math import fsum, isclose, isfinite
from pathlib import Path
import sys
from types import SimpleNamespace

VENDOR=Path(__file__).resolve().parents[1]/'_vendor'


class _RejectCoolProp(importlib.abc.MetaPathFinder):
    def __init__(self):
        self.attempts=[]

    def find_spec(self,fullname,path=None,target=None):
        if fullname.casefold().startswith('coolprop'):
            self.attempts.append(fullname)
            raise RuntimeError('CoolProp import is forbidden in Aspen-native evaluation: '+fullname)


@contextmanager
def _native_runtime():
    loaded=[x for x in sys.modules if x.casefold().startswith('coolprop')]
    if loaded:
        raise RuntimeError('Use a clean Python process without CoolProp loaded.')
    guard=_RejectCoolProp();sys.meta_path.insert(0,guard)
    path=str(VENDOR)
    inserted=path not in sys.path
    if inserted:sys.path.insert(0,path)
    try:
        import tespy
        if not Path(tespy.__file__).resolve().is_relative_to(VENDOR.resolve()):
            raise RuntimeError('An unverified external TESPy module was already imported; start a clean process.')
        from tespy.components import HeatBus,HeatExchanger,SimpleHeatExchanger,Compressor,Turbine,Pump,Merge,DropletSeparator
        yield dict(tespy=tespy,HeatBus=HeatBus,HeatExchanger=HeatExchanger,SimpleHeatExchanger=SimpleHeatExchanger,
                   Compressor=Compressor,Turbine=Turbine,Pump=Pump,Merge=Merge,DropletSeparator=DropletSeparator),guard
    finally:
        sys.meta_path.remove(guard)
        if inserted:sys.path.remove(path)


def _finite(value,label):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not isfinite(value):
        raise ValueError(f'{label}: finite imported Aspen value required')
    return float(value)


def _port(rate_kw):
    # Native equation expects an SI rate container. This is a data adapter,
    # not a replacement for any TESPy computational method.
    return SimpleNamespace(E=SimpleNamespace(val_SI=rate_kw*1000.))


def _material_port(s):
    m=_finite(s.get('mass_kg_h'),s['id']+' mass')/3600.
    H=_finite(s.get('enthalpy_kw'),s['id']+' enthalpy')*1000.
    if m<=0:
        raise ValueError('Native specific-enthalpy component check requires positive mass flow.')
    return SimpleNamespace(m=SimpleNamespace(val_SI=m),h=SimpleNamespace(val_SI=H/m),
                           label=s['id'])


def _rows(value):
    return list(value.values()) if isinstance(value,dict) else list(value or [])


def _bus_residual(native,name,incoming,outgoing):
    bus=native['HeatBus']('Aspen energy ledger '+name)
    bus.heat_inl=[_port(x) for x in incoming]
    bus.heat_outl=[_port(x) for x in outgoing]
    return float(bus.energy_balance_func())/1000.


def _component_check(native,block,raw_streams,transfer_rates):
    """Native component equations only for compatible imported port structures."""
    ins=[raw_streams[x] for x in block['inlets']]
    outs=[raw_streams[x] for x in block['outlets']]
    typ=block.get('type','').upper()
    if any(not isinstance(s.get('mass_kg_h'),(int,float)) or not isfinite(s['mass_kg_h']) for s in ins+outs):
        return {'status':'not_applicable','reason':'Native equipment check needs finite material mass flows.'}
    kind=None;note='No compatible native equipment equation selected; native boundary ledger remains available.'
    ordered=outs
    if typ in {'HEATER','COMPR','PUMP'} and len(ins)==len(outs)==1:
        if not isclose(ins[0]['mass_kg_h'],outs[0]['mass_kg_h'],rel_tol=1e-9,abs_tol=1e-7):
            return {'status':'not_applicable','reason':'Single-stream native component assumes equal inlet/outlet mass flow.'}
        if typ=='HEATER':kind='SimpleHeatExchanger'
        elif typ=='PUMP':kind='Pump'
        else:kind='Compressor' if fsum(transfer_rates)>=0 else 'Turbine'
    elif typ=='HEATX' and len(ins)==len(outs)==2 and not any(abs(x)>1e-9 for x in transfer_rates):
        # Match conserved channel mass rates; ambiguity is harmless only when
        # both inlet flows are equal, because the native equation is symmetric.
        for candidate in (outs,list(reversed(outs))):
            if all(isclose(a['mass_kg_h'],b['mass_kg_h'],rel_tol=1e-9,abs_tol=1e-7) for a,b in zip(ins,candidate)):
                ordered=candidate;kind='HeatExchanger';break
    elif typ=='MIXER' and len(outs)==1 and not any(abs(x)>1e-9 for x in transfer_rates):
        kind='Merge'
    elif typ=='FLASH2' and len(ins)==1 and len(outs)==2 and not any(abs(x)>1e-9 for x in transfer_rates):
        kind='DropletSeparator'
    if kind is None:return {'status':'not_applicable','reason':note}
    try:
        obj=native[kind]('Aspen component '+block['id'])
        obj.inl=[_material_port(s) for s in ins]
        obj.outl=[_material_port(s) for s in ordered]
        applied=fsum(transfer_rates)*1000.
        if kind=='SimpleHeatExchanger':obj.Q.val_SI=applied
        elif kind in {'Compressor','Turbine','Pump'}:obj.P.val_SI=applied
        residual=float(obj.energy_balance_func())/1000.
        # TESPy machine/HX convention is out-minus-in; normalize for comparison
        # to the Aspen control-volume convention in-minus-out.
        normalized=residual if kind in {'Merge','DropletSeparator'} else -residual
        return {'status':'evaluated_imported_states','native_class':kind,
                'native_residual_kw':residual,'energy_residual_kw':normalized,
                'method':type(obj).energy_balance_func.__module__+'.'+type(obj).energy_balance_func.__qualname__,
                'boundary_note':'For machines, P is net transferred energy to fluid after the documented drive heat loss; no efficiency is refitted.'}
    except ValueError as exc:
        return {'status':'not_applicable','reason':str(exc)}


def run_tespy(model,stream_results,balance_result):
    """Return all imported-state energy residuals computed by native TESPy.

    Boundary membership and crossing-transfer metadata come from the project
    topology adapter. No previously calculated residual is an input to TESPy.
    Material enthalpy rates and individual Q/W values are taken from Aspen data.
    """
    raw={s['id']:s for s in model['streams']}
    raw_blocks={b['id']:b for b in model['blocks']}
    block_rows=_rows(balance_result['blocks'])
    group_rows=_rows(balance_result.get('groups',[]))
    rows=block_rows+group_rows+[balance_result['plant']]
    results={};component_checks={}
    with _native_runtime() as (native,guard):
        for boundary in rows:
            name=boundary['id']
            if name in results:raise ValueError('Duplicate TESPy boundary '+name)
            block=raw_blocks.get(name)
            try:
                key='transfers' if block else 'crossing_transfers'
                source=block if block else boundary
                if key not in source or not isinstance(source[key],list):
                    raise ValueError(name+': explicit transfer inventory is required; declare [] only for an adiabatic/work-free boundary.')
                transfers=source[key]
                incoming=[_finite(raw[s].get('enthalpy_kw'),s+' enthalpy') for s in boundary['inlets']]
                outgoing=[_finite(raw[s].get('enthalpy_kw'),s+' enthalpy') for s in boundary['outlets']]
                rates=[_finite(t.get('rate_kw'),name+'/'+t.get('id','?')+' transfer') for t in transfers]
                incoming += [x for x in rates if x>=0]
                outgoing += [-x for x in rates if x<0]
                residual=_bus_residual(native,name,incoming,outgoing)
                results[name]={'energy_residual_kw':residual,'status':'evaluated_imported_states',
                    'method':'tespy.components.energy._bus._EnergyBus.energy_balance_func',
                    'input_port_count':len(incoming),'output_port_count':len(outgoing),
                    'boundary_kind':'block' if block else ('plant' if name=='plant' else 'package'),
                    'enthalpy_scope':'Whole material stream including all supplied Aspen substreams',
                    'interpretation':'Generic control-volume energy ledger; not a solved TESPy equipment model.'}
                if block:component_checks[name]=_component_check(native,block,raw,rates)
            except (ValueError,KeyError) as exc:
                results[name]={'energy_residual_kw':None,'status':'unavailable','reason':str(exc)}
        manifest=json.loads((VENDOR/'tespy/VENDOR_MANIFEST.json').read_text(encoding='utf8'))
        native_version=native['tespy'].__version__
    return {'engine':'TESPy','version':native_version,'boundaries':results,
        'component_checks':component_checks,
        'coverage':{'boundaries':len(results),'evaluated':sum(v['energy_residual_kw'] is not None for v in results.values()),
                    'native_equipment_checks':sum(v['status']=='evaluated_imported_states' for v in component_checks.values())},
        'provenance':{'source_url':manifest['source_url'],'upstream_wheel_sha256':manifest['sha256'],
            'vendor_profile':manifest['profile'],'patched_import_files':manifest['changed_files'],
            'equations_byte_identical':all(x['byte_identical'] for x in manifest['unchanged_native_equation_files']),
            'CoolProp_import_attempts':guard.attempts,'CoolProp_loaded':any(x.casefold().startswith('coolprop') for x in sys.modules),
            'process_network_solved':False,'physical_properties_recalculated':False,
            'exergy_role':'TESPy 0.12 delegates exergy analysis to ExerPy; this adapter evaluates energy residuals.'}}
