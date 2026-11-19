"""Independent, composition-preserving estimates of physical exergy.

Peng–Robinson uses ChEDL's fluid constants and heat capacities, classical mixing
with zero binary interaction parameters, and at most one vapor and one liquid.
Both end states of each difference use the same property model.
"""
from functools import lru_cache
from math import fsum, isclose, isfinite


class PropertyUnavailable(ValueError):
    """No defensible estimate could be made with the selected property model."""


@lru_cache(maxsize=512)
def _identify(identifier):
    from chemicals.elements import simple_formula_parser
    from chemicals.identifiers import search_chemical

    for spelling in dict.fromkeys([identifier, identifier.replace("-", " ")]):
        try:
            record = search_chemical(spelling)
            return record.CASs, simple_formula_parser(record.formula)
        except ValueError:
            continue
    raise PropertyUnavailable(f"Unrecognized chemical identifier {identifier!r}; supply the component CAS number.")


@lru_cache(maxsize=64)
def _flash_model(cas_numbers):
    from thermo import CEOSGas, CEOSLiquid, ChemicalConstantsPackage, FlashPureVLS, FlashVL, PRMIX

    constants, properties = ChemicalConstantsPackage.from_IDs(cas_numbers)
    eos = dict(Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas)
    if any(x is None for values in eos.values() for x in values):
        raise PropertyUnavailable("Peng–Robinson critical properties or acentric factor missing.")
    gases = properties.HeatCapacityGases
    if any(cp.method is None for cp in gases):
        raise PropertyUnavailable("Missing ideal-gas heat capacity.")
    gas = CEOSGas(PRMIX, eos_kwargs=eos, HeatCapacityGases=gases)
    liquid = CEOSLiquid(PRMIX, eos_kwargs=eos, HeatCapacityGases=gases)
    if len(cas_numbers) == 1:
        flash = FlashPureVLS(constants, properties, gas=gas, liquids=[liquid], solids=[])
    else:
        flash = FlashVL(constants, properties, gas=gas, liquid=liquid)
        flash.PT_SS_MAXITER = 200
        flash.PT_SS_POLISH_MAXITER = 100
    return flash, constants, gases


def is_graphite(component):
    return (component.get("atoms") == {"C": 1}
            and (component.get("cas") == "7782-42-5"
                 or component.get("name", "").casefold() in {"carbon-graphite", "graphite"}))


@lru_cache(maxsize=1)
def _graphite_cp():
    from thermo import HeatCapacitySolid
    return HeatCapacitySolid(CASRN="7782-42-5", method="JANAF")


@lru_cache(maxsize=256)
def _water_specific(t, p, vapor_fraction):
    from iapws import IAPWS97
    pressure_mpa = p * 1.e-6
    if vapor_fraction is not None and 0.000611657 <= pressure_mpa < 22.064:
        saturated = IAPWS97(P=pressure_mpa, x=vapor_fraction)
        if abs(t - saturated.T) <= 0.4:
            return saturated.h, saturated.s
        if 0 < vapor_fraction < 1 and min(vapor_fraction, 1 - vapor_fraction) > 1e-9:
            raise PropertyUnavailable("Water quality and temperature disagree with the IF97 saturation curve.")
    state = IAPWS97(P=pressure_mpa, T=t)
    return state.h, state.s


@lru_cache(maxsize=512)
def _molar_availability(ids, zs, t, p, t0, vapor_fraction=None):
    flash, _, _ = _flash_model(ids)
    kwargs = {"P": p, "T": t}
    if len(zs) == 1 and vapor_fraction is not None and 1e-9 < vapor_fraction < 1 - 1e-9:
        kwargs = {"P": p, "VF": vapor_fraction}
    state = flash.flash(zs=list(zs), **kwargs)
    if abs(state.T - t) > 2.0:
        raise PropertyUnavailable("Pure-fluid saturation estimate disagrees with the operating temperature.")
    if not all(0 <= b <= 1 for b in state.betas):
        raise PropertyUnavailable("Flash produced an invalid phase fraction.")
    if len(state.phases) > 1:
        errors = [abs(fsum(b * phase.zs[i] for b, phase in zip(state.betas, state.phases)) - z)
                  for i, z in enumerate(zs)]
        if max(errors) > 2e-7:
            raise PropertyUnavailable("Estimated flash composition fails material closure.")
    value = state.H() - t0 * state.S()
    if not isfinite(value):
        raise PropertyUnavailable("Nonfinite flash result.")
    return value


def _wilson_availability(flash, constants, zs, t, p, t0):
    """Use Wilson K-values with a material-balanced split, then evaluate PR H/S."""
    from chemicals.flash_basic import Wilson_K_value
    from chemicals.rachford_rice import flash_inner_loop

    ratios = [Wilson_K_value(t, p, tc, pc, omega) for tc, pc, omega in
              zip(constants.Tcs, constants.Pcs, constants.omegas)]
    if fsum(z * k for z, k in zip(zs, ratios)) <= 1:
        phases = [(1., flash.liquids[0], zs)]
    elif fsum(z / k for z, k in zip(zs, ratios)) <= 1:
        phases = [(1., flash.gas, zs)]
    else:
        vf, xs, ys = flash_inner_loop(zs, ratios)
        if not 0 <= vf <= 1 or any(abs(vf * y + (1-vf) * x - z) > 2e-7 for z, x, y in zip(zs, xs, ys)):
            raise PropertyUnavailable("Wilson phase split failed composition closure.")
        phases = [(vf, flash.gas, ys), (1 - vf, flash.liquids[0], xs)]
    values = []
    for amount, phase, composition in phases:
        state = phase.to(T=t, P=p, zs=composition)
        values.append(amount * (state.H() - t0 * state.S()))
    return fsum(values)


@lru_cache(maxsize=256)
def _fluid_delta(ids, zs, t, p, t0, p0, vapor_fraction):
    """Cache only molar differences, with both state specifications in the key."""
    notes = []
    flash, constants, cp_gases = _flash_model(ids)
    try:
        operating = _molar_availability(ids, zs, t, p, t0, vapor_fraction)
        reference = _molar_availability(ids, zs, t0, p0, t0, None)
    except (MemoryError, RecursionError):
        raise
    except Exception as error:
        if len(zs) == 1:
            raise PropertyUnavailable(f"Pure-fluid flash failed: {error}") from error
        operating = _wilson_availability(flash, constants, list(zs), t, p, t0)
        reference = _wilson_availability(flash, constants, list(zs), t0, p0, t0)
        if not all(isfinite(v) for v in (operating, reference)):
            raise PropertyUnavailable("Wilson phase approximation produced nonfinite availability.")
        method = "thermo_pr_wilson_k"
        notes.append(f"PR flash failed ({type(error).__name__}). Wilson K-values approximate phase partition at both states.")
    else:
        method = "thermo_peng_robinson"
    notes.append("One vapor and one liquid phase allowed, zero binary interaction parameters.")
    if any(z > 1e-7 and tm is not None and min(t0, t) < tm for z, tm in zip(zs, constants.Tms)):
        notes.append("Solid precipitation is omitted. A constituent's pure melting point exceeds a flash temperature.")
    if any(z > 1e-7 and (cp.Tmin is not None and min(t0, t) < cp.Tmin or
                          cp.Tmax is not None and max(t0, t) > cp.Tmax) for z, cp in zip(zs, cp_gases)):
        notes.append("At least one ideal-gas heat-capacity correlation is extrapolated.")
    return operating - reference, method, tuple(notes)


def _estimate(stream, active, t0, p0, formation):
    t, p = stream.get("temperature_k"), stream.get("pressure_pa")
    if t is None or p is None:
        raise PropertyUnavailable("Operating temperature and pressure are required for a model estimate.")
    if t <= 0 or p <= 0:
        raise PropertyUnavailable("Operating temperature and pressure must be positive.")
    notes = []
    solids = fluids = 0.0
    rates = {}
    for component, n in active:
        if is_graphite(component):
            native = (stream.get("substreams") or {}).get("CISOLID", {})
            graphite = native.get("component_flows_kmol_h", {})
            if (formation and set(graphite) == {component["id"]}
                    and isclose(graphite[component["id"]], n, rel_tol=1e-11, abs_tol=0.)
                    and all(native.get(k) is not None for k in ("enthalpy_kw", "entropy_kw_k"))):
                solids += native["enthalpy_kw"] - t0 * native["entropy_kw_k"]
                notes.append("Graphite uses the saved Aspen solid substream and elemental reference.")
                continue
            cp = _graphite_cp()
            solids += n / 3600 * (cp.T_dependent_property_integral(t0, t)
                                  - t0 * cp.T_dependent_property_integral_over_T(t0, t))
            notes.append("Graphite uses JANAF solid heat capacity with pressure dependence neglected.")
            continue
        if component.get("type", "").lower() not in {"", "conventional", "fluid"}:
            raise PropertyUnavailable(f"Unsupported component type for {component['id']}.")
        identifier = component.get("cas") or component.get("name") or component["id"]
        cas, atoms = _identify(identifier.lower())
        if component.get("atoms") and atoms != component["atoms"]:
            raise PropertyUnavailable(f"Chemical identifier and elemental formula disagree for {component['id']}.")
        rates[cas] = rates.get(cas, 0) + n
    if not rates:
        return solids, "janaf_graphite", notes

    has_graphite = any(is_graphite(c) for c, _ in active)
    vf = None if has_graphite else stream.get("vapor_fraction")
    if vf is not None:
        vf = min(1., max(0., vf))
    if set(rates) == {"7732-18-5"}:
        h0, s0 = _water_specific(t0, p0, None)
        h, s = _water_specific(t, p, vf)
        delta = (h - h0) - t0 * (s - s0)
        # IF97 returns specific kJ/kg and kJ/(kg K). Use saved mass for pure water.
        mass = stream.get("mass_kg_h") if len(active) == 1 else None
        if mass is None:
            mass = rates["7732-18-5"] * 18.015268
        fluids = mass / 3600 * delta
        method = "iapws_if97"
    else:
        ids = tuple(sorted(rates))
        flow = fsum(rates.values())
        zs = tuple(rates[c] / flow for c in ids)
        delta, method, model_notes = _fluid_delta(ids, zs, t, p, t0, p0, vf if len(zs) == 1 else None)
        fluids = delta * flow / 3600
        notes.extend(model_notes)
    if has_graphite:
        method += "_with_graphite"
    return fsum([solids, fluids]), method, notes


def estimate_physical(stream, active, t0, p0, *, formation=False):
    """Return a signed kW estimate, property-model ID, and model limitations."""
    from chemicals.exceptions import PhaseCountReducedError, PhaseExistenceImpossible, TrivialSolutionError
    from fluids.numerics import NoSolutionError, NotBoundedError, OscillationError, UnconvergedError

    try:
        value, method, notes = _estimate(stream, active, t0, p0, formation)
        if not isfinite(value):
            raise PropertyUnavailable("Property evaluation returned a nonfinite physical exergy.")
        return value, method, notes
    except PropertyUnavailable:
        raise
    except (ArithmeticError, ValueError, TypeError, RuntimeError, KeyError, PhaseCountReducedError,
            PhaseExistenceImpossible, TrivialSolutionError,
            NoSolutionError, NotBoundedError, OscillationError, UnconvergedError) as exc:
        raise PropertyUnavailable(f"Property evaluation failed ({type(exc).__name__}): {exc}") from exc
