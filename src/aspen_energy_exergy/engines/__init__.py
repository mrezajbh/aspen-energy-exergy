"""Explicit imported-state evaluation by pinned ExerPy and TESPy source releases.

These backends evaluate balances, not an Aspen process resimulation or an
independent fluid-property prediction. See the vendored source provenance.
"""
from contextlib import contextmanager
import importlib.abc
from math import fsum, isclose
import sys


@contextmanager
def no_coolprop():
    """Fail on an accidental CoolProp import; do not provide a replacement module."""
    existing=[n for n in sys.modules if n == 'CoolProp' or n.startswith('CoolProp.')]
    if existing:
        raise RuntimeError('This execution requires a fresh Python process without CoolProp loaded.')
    attempts=[]
    class Guard(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path=None, target=None):
            if fullname == 'CoolProp' or fullname.startswith('CoolProp.'):
                attempts.append(fullname)
                raise ImportError('CoolProp is prohibited in the Aspen imported-state calculation.')
    guard=Guard();sys.meta_path.insert(0,guard)
    try:
        yield attempts
        if attempts:
            raise RuntimeError('A backend attempted to import CoolProp: '+', '.join(attempts))
    finally:
        sys.meta_path.remove(guard)


def run_external_engines(model, stream_results, balance_result):
    """Execute both upstream engines and compare identical boundary definitions.

No result is called an independent property validation. Missing states remain
unavailable. Negative signed destruction and residuals are retained.
"""
    with no_coolprop() as attempts:
        from .exerpy_engine import run_exerpy
        from .tespy_engine import run_tespy
        ex=run_exerpy(model,stream_results,balance_result)
        en=run_tespy(model,stream_results,balance_result)
    original={b['id']:b for b in balance_result['blocks']}
    original.update(balance_result.get('groups',{}))
    original['plant']=balance_result['plant']
    comparisons=[]
    for engine,report,field in [('ExerPy',ex,'exergy_destruction_kw'),('TESPy',en,'energy_residual_kw')]:
        for bid,record in report['boundaries'].items():
            expected=original[bid].get(field)
            actual=record.get(field)
            delta=None if expected is None or actual is None else actual-expected
            matched=None if delta is None else isclose(actual,expected,rel_tol=1e-9,abs_tol=1e-6)
            comparisons.append({'boundary':bid,'engine':engine,'quantity':field,
                'previous_formula_kw':expected,'library_kw':actual,'difference_kw':delta,
                'agreement_within_tolerance':matched,'method':record.get('method'),
                'status':record.get('status')})
    return {'schema':'aspen-external-balances/v1','ExerPy':ex,'TESPy':en,
        'comparisons':comparisons,
        'coolprop':{'import_attempts':attempts,'loaded':any(n=='CoolProp' or n.startswith('CoolProp.') for n in sys.modules)},
        'scope':'Imported Aspen states and declared boundary transfers; upstream balance evaluation, not independent property validation or a network solve.',
        'tolerance':{'absolute_kw':1e-6,'relative':1e-9}}


def promote_external_results(report):
    """Publish verified upstream balances in place and return ``report``.

    ``report`` must contain ``streams``, ``blocks``, ``groups``, ``plant`` and
    the result of :func:`run_external_engines` in ``external_engines``. All
    comparisons are checked before mutation. Missing upstream results remain
    null, rather than silently falling back to a built-in value. Supporting
    material/transfer inventories still come from the supplied Aspen states.
    This is imported-state balance evaluation, not thermodynamic resimulation.
    """
    external = report['external_engines']
    failed = [r for r in external['comparisons'] if r['agreement_within_tolerance'] is False]
    if failed:
        raise ValueError('External balance disagreement: ' + str(failed))
    blocks = {r['id']: r for r in report['blocks']}
    boundaries = {**blocks, **report.get('groups', {}), 'plant': report['plant']}
    for engine in ('ExerPy', 'TESPy'):
        if set(external[engine]['boundaries']) != set(boundaries):
            raise ValueError(f'{engine} boundary inventory differs from the report')
    if set(external['ExerPy']['stream_results']) != set(report['streams']):
        raise ValueError('ExerPy stream inventory differs from the report')

    for bid, target in boundaries.items():
        for engine, field in (('ExerPy', 'exergy_destruction_kw'), ('TESPy', 'energy_residual_kw')):
            result = external[engine]['boundaries'][bid]
            target[field] = result[field]
            target[field + '_engine'] = engine
            target[field + '_provenance'] = {
                'engine': engine, 'version': external[engine]['version'],
                'method': result['method'], 'status': result['status'],
                'basis': 'Imported Aspen states and declared boundary transfers',
            }
        destruction, residual = target['exergy_destruction_kw'], target['energy_residual_kw']
        target['exergy_balance_status'] = ('incomplete' if destruction is None else
                                          'negative' if destruction < -1e-3 else 'available')
        tolerance = target.get('energy_balance_tolerance_kw')
        target['energy_balance_status'] = ('incomplete' if residual is None or tolerance is None else
                                          'closed' if abs(residual) <= tolerance else 'residual')
        # These are identities/diagnostics of the published balances, not new
        # upstream property predictions or substitutes for upstream equations.
        lost_work = target.get('t0_entropy_generation_kw')
        reference = target.get('element_reference_residual_kw')
        target['exergy_minus_t0_sgen_kw'] = (None if destruction is None or lost_work is None
                                            else destruction - lost_work)
        target['exergy_identity_residual_kw'] = (None if None in (destruction, lost_work, residual, reference)
                                                else destruction - fsum((lost_work, residual, reference)))

    for sid, result in external['ExerPy']['stream_results'].items():
        stream = report['streams'][sid]
        stream['total_exergy_kw'] = result['total_exergy_kw']
        stream['total_exergy_aggregation'] = result['method']
        stream['total_exergy_kw_provenance'] = {
            'engine': 'ExerPy' if result['status'] == 'reaggregated_from_supplied_split' else 'imported',
            'method': result['method'], 'status': result['status'],
            'properties': 'Supplied Aspen states and established elemental-reference reconstruction',
        }

    def complete_sum(values):
        values = list(values)
        return None if None in values else fsum(values)

    for target in [*report.get('groups', {}).values(), report['plant']]:
        destruction = complete_sum(blocks[bid]['exergy_destruction_kw'] for bid in target['members'])
        residual = complete_sum(blocks[bid]['energy_residual_kw'] for bid in target['members'])
        target['block_exergy_destruction_sum_kw'] = destruction
        combined = complete_sum((destruction, target.get('internal_link_exergy_destruction_kw')))
        target['exergy_additivity_residual_kw'] = (None if None in (target['exergy_destruction_kw'], combined)
                                                  else target['exergy_destruction_kw'] - combined)
        target['energy_additivity_residual_kw'] = (None if None in (target['energy_residual_kw'], residual)
                                                  else target['energy_residual_kw'] - residual)
    plant = report['plant']
    if plant.get('product_streams'):
        plant['product_exergy_kw'] = complete_sum(report['streams'][sid]['total_exergy_kw']
                                                   for sid in plant['product_streams'])
        incoming = plant.get('exergy_input_kw')
        plant['product_exergy_fraction_of_input'] = (None if plant['product_exergy_kw'] is None or
                                                    incoming is None or incoming <= 0 else
                                                    plant['product_exergy_kw'] / incoming)
    report['balance_engine'] = 'exerpy-tespy'
    report.setdefault('software_versions', {}).update({engine: external[engine]['version']
                                                       for engine in ('ExerPy', 'TESPy')})
    return report
