"""ExerPy 0.2.2 balances of Aspen-exported states, without a property backend.

The upstream numerical routines are unmodified. A documented vendor patch only
moves optional CoolProp imports into unused property-conversion functions.
Stream states and chemical-exergy references are inputs, not predictions made by
ExerPy. A control-volume balance is distinct from a native equipment model.
"""
from __future__ import annotations

from copy import deepcopy
from contextlib import contextmanager
import importlib
import logging
from math import isclose, isfinite
from pathlib import Path
import sys


def _load():
    vendor = Path(__file__).resolve().parents[1] / "_vendor"
    if str(vendor) not in sys.path:
        sys.path.insert(0, str(vendor))
    package = importlib.import_module("exerpy")
    if Path(package.__file__).resolve().parent != vendor / "exerpy":
        raise RuntimeError("Expected the audited ExerPy vendor build, not an unverified installation")
    if any(n.lower() == "coolprop" or n.lower().startswith("coolprop.") for n in sys.modules):
        raise RuntimeError("CoolProp must not be loaded for the Aspen-only analysis")
    return package, importlib.import_module("exerpy.functions")


@contextmanager
def _messages():
    records = []
    class Capture(logging.Handler):
        def emit(self, record):
            if record.levelno >= logging.WARNING:
                records.append(record.getMessage())
    logger = logging.getLogger("exerpy")
    handler = Capture(); previous = logger.propagate
    logger.addHandler(handler); logger.propagate = False
    try:
        yield records
    finally:
        logger.removeHandler(handler); logger.propagate = previous


def _value(value):
    return float(value) if value is not None and isfinite(value) else None


def _material(stream, result):
    """Translate rates in kW/kg h to ExerPy's W, kg/s and J/kg."""
    mass = stream.get("mass_kg_h")
    if mass is None or mass < 0:
        raise ValueError(f"{stream['id']}: valid mass flow required")
    m = mass / 3600
    row = dict(name=stream["id"], kind="material", m=m,
        T=stream.get("temperature_k"), p=stream.get("pressure_pa"),
        source_component=stream.get("from_block"), target_component=stream.get("to_block"),
        source_connector=0, target_connector=0)
    for field, key in [("enthalpy_kw", "h"), ("entropy_kw_k", "s"),
                       ("physical_exergy_kw", "e_PH"), ("chemical_exergy_kw", "e_CH")]:
        rate = result.get(field, stream.get(field))
        row[key] = rate * 1000 / m if m and rate is not None else (0.0 if not m else None)
    if result.get("total_exergy_kw") is not None:
        row["E"] = result["total_exergy_kw"] * 1000
    return row


def _boundary(package, name, scope, boundary, transfers, materials, t0, p0):
    connections = {}; inlet_names = []; outlet_names = []
    missing = []
    for side, ids, names in [("in", boundary["inlets"], inlet_names), ("out", boundary["outlets"], outlet_names)]:
        for sid in ids:
            conn = deepcopy(materials[sid]); key = f"material:{side}:{sid}"
            if conn.get("E") is None: missing.append(sid)
            conn.update(source_component=None if side == "in" else name,
                        target_component=name if side == "in" else None)
            connections[key] = conn; names.append(key)
    for n, transfer in enumerate(transfers):
        exergy = transfer.get("exergy_kw")
        if exergy is None:
            missing.append(f"{transfer['block']}/{transfer['id']} heat/work exergy")
            continue
        key = f"transfer:{n}:{transfer['block']}:{transfer['id']}"
        positive = exergy >= 0
        connections[key] = dict(name=key, kind="power" if transfer["kind"] == "work" else "heat",
            E=abs(exergy) * 1000, energy_flow=abs(transfer["rate_kw"]) * 1000,
            source_component=None if positive else name, target_component=name if positive else None,
            source_connector=0, target_connector=0,
            imported_heat_exergy_basis=transfer.get("basis"))
        (inlet_names if positive else outlet_names).append(key)
    method = "ExerPy.ExergyAnalysis.analyse; external control-volume ledger"
    base = dict(id=name, scope=scope, method=method, internal_equipment_model=False,
                exergy_destruction_kw=None, exergy_input_kw=None, exergy_output_kw=None,
                reference_exergy_destruction_kw=boundary.get("exergy_destruction_kw"))
    if missing:
        return {**base, "status":"missing_input", "missing":missing}
    # An empty internal component collection is intentional: this is a control
    # volume with externally supplied material/heat/work connections. No fake
    # reactor, separator, or black-box equipment class is registered.
    with _messages() as warnings:
        analysis = package.ExergyAnalysis({}, connections, t0, p0, split_physical_exergy=False)
        analysis.analyse({"inputs":inlet_names}, {"inputs":outlet_names})
    destroyed = _value(analysis.E_D)
    previous = boundary.get("exergy_destruction_kw")
    return {**base, "status":"calculated", "exergy_destruction_kw":destroyed / 1000,
            "exergy_input_kw":analysis.E_F / 1000, "exergy_output_kw":analysis.E_P / 1000,
            "difference_kw":None if previous is None else destroyed / 1000 - previous,
            "connection_count":len(connections), "warnings":warnings,
            "warning_context":"No internal component ledger is constructed for this control-volume audit; the component-sum warning is expected. Negative-destruction warnings remain substantive.",
            "heat_exergy_source":"Imported signed heat exergy from documented heat-boundary temperatures; no ExerPy fluid-difference substitution"}


def _same_flow_recipe(a, b):
    if not isclose(a.get("mass_kg_h",0),b.get("mass_kg_h",0),rel_tol=1e-6,abs_tol=1e-5):return False
    fa,fb=a.get("component_flows_kmol_h",{}),b.get("component_flows_kmol_h",{})
    return all(isclose(fa.get(k,0),fb.get(k,0),rel_tol=1e-6,abs_tol=1e-7) for k in fa.keys()|fb.keys())


def _component(package, block, boundary, transfers, materials, streams, t0, p0):
    bid=block['id']; kind=block.get('type'); native=None
    ins=list(block['inlets']); outs=list(block['outlets'])
    if any(t.get('rate_kw') is None for t in transfers if t['kind']=='work'):
        return dict(id=bid,aspen_type=kind,native_class=None,status='missing_work_input',
                    comparable_to_full_boundary=False,exergy_destruction_kw=None,
                    reason='Known shaft-work values are required for the native class')
    work=sum(t['rate_kw'] for t in transfers if t['kind']=='work')
    if kind=='COMPR':native='Compressor' if work>=0 else 'Turbine'
    elif kind=='PUMP':native='Pump'
    elif kind=='VALVE':native='Valve'
    elif kind=='MIXER':native='Mixer'
    elif kind in ('FSPLIT','SSPLIT'):native='Splitter'
    elif kind=='FLASH2':native='FlashTank'
    elif kind=='HEATX':native='HeatExchanger'
    result=dict(id=bid, aspen_type=kind, native_class=native, status='unsupported_native_equipment',
                comparable_to_full_boundary=False, exergy_destruction_kw=None,
                full_boundary_exergy_destruction_kw=boundary.get('exergy_destruction_kw'))
    if native=='Mixer':
        # A zero-flow uninitialized Aspen inlet has no temperature. It contributes
        # exactly zero to the complete ledger and is omitted only from the native
        # mixing-temperature classification, without inventing a state for it.
        inactive=[s for s in ins if materials[s]['m']==0 and materials[s].get('T') is None]
        ins=[s for s in ins if s not in inactive]
        result['omitted_zero_flow_native_ports']=inactive
    if native is None:
        result['reason']='No applicable built-in ExerPy equipment class; full control-volume ledger is still evaluated.'
        return result
    if any(materials[s].get('e_PH') is None or materials[s].get('e_CH') is None for s in ins+outs):
        return {**result,'status':'missing_physical_chemical_split','reason':'Current Aspen/reference-matched split unavailable; total exergy is not relabelled as physical exergy.'}
    if native in ('Compressor','Pump','Valve') and (len(ins)!=1 or len(outs)!=1):
        return {**result,'reason':'Native class requires one material inlet and outlet'}
    if native=='HeatExchanger':
        if len(ins)!=2 or len(outs)!=2:return {**result,'reason':'Native two-sided class needs two inlet/outlet streams'}
        ins.sort(key=lambda s:materials[s]['T'],reverse=True)
        pairs=[[o for o in outs if _same_flow_recipe(streams[i],streams[o])] for i in ins]
        if any(len(p)!=1 for p in pairs) or pairs[0]==pairs[1]:
            return {**result,'reason':'Heat-exchanger side pairing is not unambiguous from mass/component flows'}
        outs=[p[0] for p in pairs]
    connections={}
    for side,ids in [('in',ins),('out',outs)]:
        for index,sid in enumerate(ids):
            row=deepcopy(materials[sid]);row.update(source_component=None if side=='in' else bid,
                target_component=bid if side=='in' else None,source_connector=index,target_connector=index)
            connections[f'{side}:{sid}']=row
    if native in ('Compressor','Pump','Turbine'):
        positive=work>=0
        connections['shaft']=dict(name='shaft',kind='power',energy_flow=abs(work)*1000,E=abs(work)*1000,
            source_component=None if positive else bid,target_component=bid if positive else None,
            source_connector=len(outs),target_connector=len(ins))
    with _messages() as warnings:
        try:
            analysis=package.ExergyAnalysis({native:{bid:{'name':bid}}},connections,t0,p0,split_physical_exergy=False)
            component=analysis.components[bid]
            component.calc_exergy_balance(t0,p0,False)
            val=_value(component.E_D)
        except (ValueError,KeyError,TypeError,AttributeError,ZeroDivisionError) as error:
            return {**result,'status':'native_class_rejected','reason':str(error),'warnings':warnings}
    reasons=[]
    if any(t.get('exergy_kw') is None for t in transfers if t['kind']=='heat'):
        reasons.append('External heat exergy is unknown; native class cannot validate the complete boundary')
    heat=sum(t.get('exergy_kw') or 0 for t in transfers if t['kind']=='heat')
    if abs(heat)>1e-6:reasons.append('Native equipment class omits nonzero external heat exergy')
    if native not in ('Compressor','Pump','Turbine') and abs(work)>1e-6:
        reasons.append('Native equipment class omits external work')
    ch=sum(materials[s]['m']*materials[s]['e_CH'] for s in ins)-sum(materials[s]['m']*materials[s]['e_CH'] for s in outs)
    if native not in ('Mixer','FlashTank') and abs(ch)>1e-3:
        reasons.append('Native equipment class excludes non-cancelling chemical exergy in supplied states')
    if abs(sum(materials[s]['m'] for s in ins)-sum(materials[s]['m'] for s in outs))>1e-7:
        reasons.append('Material mass imbalance changes native-class assumptions')
    return {**result,'status':'calculated' if val is not None else 'native_result_unavailable',
            'exergy_destruction_kw':None if val is None else val/1000,
            'difference_kw':None if val is None or boundary.get('exergy_destruction_kw') is None else val/1000-boundary['exergy_destruction_kw'],
            'comparable_to_full_boundary':not reasons,'limitations':reasons,'warnings':warnings,
            'native_fuel_exergy_kw':None if _value(component.E_F) is None else component.E_F/1000,
            'native_product_exergy_kw':None if _value(component.E_P) is None else component.E_P/1000,
            'native_efficiency':_value(component.epsilon)}


def run_exerpy(model, stream_results, balance_result):
    """Run audited ExerPy stream aggregation, boundary balances and native checks.

    ``stream_results`` is keyed by stream ID; ``balance_result`` is the existing
    balance dictionary carrying block/group ports and signed heat/work legs.
    ExerPy recomputes all sums. Prior destruction values are used only to report
    differences, never as inputs to the ExerPy result.
    """
    package,functions=_load()
    streams={s['id']:s for s in model['streams']}
    materials={sid:_material(stream,stream_results[sid]) for sid,stream in streams.items()}
    known={sid:row for sid,row in materials.items() if row['e_PH'] is not None and row['e_CH'] is not None}
    functions.add_total_exergy_flow({'components':{},'connections':known},False)
    results={}
    for sid,row in materials.items():
        total=row.get('E');previous=stream_results[sid].get('total_exergy_kw')
        results[sid]=dict(total_exergy_kw=None if total is None else total/1000,
            physical_exergy_kw=stream_results[sid].get('physical_exergy_kw'),
            chemical_exergy_kw=stream_results[sid].get('chemical_exergy_kw'),
            method='ExerPy.add_total_exergy_flow' if sid in known else 'Imported Aspen-based total; split unavailable',
            status='reaggregated_from_supplied_split' if sid in known else 'imported_total_no_split',
            difference_kw=None if total is None or previous is None else total/1000-previous)
    t0=298.15;p0=101325.0
    boundaries={}; checks={}; transfer_rows=balance_result['transfers']
    blocks={b['id']:b for b in model['blocks']}
    for row in balance_result['blocks']:
        legs=[t for t in transfer_rows if t['block']==row['id']]
        boundaries[row['id']]=_boundary(package,row['id'],'block',row,legs,materials,t0,p0)
        checks[row['id']]=_component(package,blocks[row['id']],row,legs,materials,streams,t0,p0)
    for name,row in balance_result.get('groups',{}).items():
        boundaries[name]=_boundary(package,name,'group',row,row['crossing_transfers'],materials,t0,p0)
    plant=balance_result['plant']
    boundaries['plant']=_boundary(package,'plant','plant',plant,plant['crossing_transfers'],materials,t0,p0)
    if any(n.lower()=='coolprop' or n.lower().startswith('coolprop.') for n in sys.modules):
        raise RuntimeError('Unexpected CoolProp import')
    return dict(engine='ExerPy',version=package.__version__,stream_results=results,boundaries=boundaries,
        component_checks=checks,provenance={
            'implementation':'Audited upstream ExerPy 0.2.2; only optional-import placement patched',
            'upstream_repository':'https://github.com/oemof/exerpy',
            'release_source':'https://pypi.org/project/exerpy/0.2.2/',
            'official_documentation':'https://exerpy.readthedocs.io/',
            'properties':'Aspen exported H/S and physical exergy; chemical exergy from existing elemental-reference reconstruction. ExerPy is not an independent property prediction.',
            'boundary_method':'ExergyAnalysis.analyse on each externally defined control volume, with no fictitious internal component models',
            'heat_exergy':'Existing documented heat-boundary temperatures and signed heat-exergy transfers supplied unchanged',
            'native_equipment_checks':'Separate from complete boundary ledger; unsupported classes and missing splits remain explicitly flagged',
            'no_coolprop_loaded':True,'temperature_k':t0,'pressure_pa':p0})
