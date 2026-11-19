"""Fuel-product efficiency on an explicitly selected lower-heating-value basis."""

from functools import lru_cache
from math import fsum, isclose, isfinite
import sys

from ._estimate import PropertyUnavailable, _identify, is_graphite
from .balances import _difference, _total


# These constituents receive no combustible-fuel credit, including moisture.
_NONFUELS = {
    "7732-18-5", "124-38-9", "7727-37-9", "7782-44-7", "7446-09-5", "7446-11-9",
    "7440-59-7", "7440-01-9", "7440-37-1", "7439-90-9", "7440-63-3",
}


def _disable_optional_coolprop():
    """Select the ChEDL-only reference-property path before constructing objects.

    thermo 0.6.x has no public backend-disable argument. Its documented source
    uses ``thermo.coolprop._has_CoolProp`` as the lazy availability cache; None
    would attempt importing the optional package while enumerating methods.
    Explicitly marking this backend unavailable prevents even that probe. No
    external module is mocked and no reference data or equations are changed.
    Fail closed if the pinned thermo-series availability contract changes.
    """
    if any(name == 'CoolProp' or name.startswith('CoolProp.') for name in sys.modules):
        raise RuntimeError('Fuel performance requires a clean process without CoolProp loaded')
    from thermo import coolprop
    if not hasattr(coolprop, '_has_CoolProp') or not callable(getattr(coolprop, 'has_CoolProp', None)):
        raise RuntimeError('Unsupported thermo optional-backend availability contract')
    coolprop._has_CoolProp = False


@lru_cache(maxsize=512)
def _standard_lhv(cas):
    from chemicals.combustion import combustion_data, combustion_stoichiometry

    if cas == "7782-42-5":
        return -combustion_data(formula="C", Hf=0).LHV, "Standard-state graphite oxidation to CO2"
    if cas in _NONFUELS:
        return 0.0, "Nonfuel constituent"
    _disable_optional_coolprop()
    from thermo import ChemicalConstantsPackage
    constants = ChemicalConstantsPackage.constants_from_IDs([cas])
    atoms = constants.atomss[0]
    if set(atoms) - {"C", "H", "O", "N", "S"}:
        raise PropertyUnavailable("Heating-value model covers CHONS fuels")
    if combustion_stoichiometry(atoms).get("O2", 0.0) >= 0:
        return 0.0, "No oxygen demand in the complete-combustion convention"
    lhv = constants.Hcs_lower[0]
    method = "ChEDL Hcs_lower, pure constituent in its standard phase"
    if lhv is None:
        # Use Hess's law when only gas formation data are available. Convert
        # the fuel to its actual standard phase before calculating combustion.
        gas_hf = constants.Hfgs[0]
        phase = constants.phase_STPs[0]
        if gas_hf is None or phase not in {"l", "s"}:
            raise PropertyUnavailable("No standard-phase heat of combustion or convertible formation data")
        _, props = ChemicalConstantsPackage.from_IDs([cas])
        props_by_phase = {"l": props.EnthalpyVaporizations, "s": props.EnthalpySublimations}
        change = props_by_phase[phase][0]
        h_change = change(298.15)
        if h_change is None or not isfinite(h_change) or h_change <= 0:
            raise PropertyUnavailable("No enthalpy of phase change at 298.15 K")
        lhv = combustion_data(formula=constants.formulas[0], Hf=gas_hf - h_change).LHV
        method = f"Estimate from gas Hf and 298.15 K phase-change enthalpy ({change.method})"
    if lhv is None or not isfinite(lhv) or lhv >= 0:
        raise PropertyUnavailable("Positive fuel lower heating value unavailable")
    return -lhv, method


def _lhv_record(component):
    cas = None
    try:
        if is_graphite(component):
            cas = "7782-42-5"
        else:
            ident = component.get("cas") or component.get("name") or component["id"]
            cas, atoms = _identify(ident.lower())
            if component.get("atoms") and atoms != component["atoms"]:
                raise PropertyUnavailable("Chemical identifier and elemental formula disagree")
        lhv, method = _standard_lhv(cas)
        return dict(cas=cas, lhv_j_mol=lhv, method=method)
    except (ArithmeticError, ValueError, RuntimeError, TypeError, KeyError) as exc:
        return dict(cas=cas, lhv_j_mol=None, method="unavailable", reason=str(exc))


def _rates(stream):
    flows = stream.get("component_flows_kmol_h")
    if flows is None:
        return {} if stream.get("molar_kmol_h") == 0 else None
    if not isinstance(flows, dict) or any(
        isinstance(v, bool) or not isinstance(v, (int, float)) or not isfinite(v) or v < 0
        for v in flows.values()
    ):
        return None
    total = stream.get("molar_kmol_h")
    if total is not None and not isclose(fsum(flows.values()), total, abs_tol=1e-10, rel_tol=1e-7):
        return None
    return flows if flows or total == 0 else None


def _selected(ids, allowed, label):
    if not isinstance(ids, (list, tuple)) or not all(isinstance(s, str) for s in ids):
        raise ValueError(f"{label} must be a sequence of external stream IDs")
    if len(set(ids)) != len(ids) or not set(ids).issubset(allowed):
        raise ValueError(f"{label} contains a duplicate or stream on the wrong boundary")
    return list(ids)


def _utility_pair(inlet, outlet, diagnostics):
    source, sink = inlet["id"], outlet["id"]
    drop = _difference(inlet.get("enthalpy_kw"), outlet.get("enthalpy_kw"))
    incoming, outgoing = _rates(inlet), _rates(outlet)
    if incoming is None or outgoing is None or drop is None:
        drop = None
        diagnostics.append(f"Missing utility flow, composition or enthalpy for {source}/{sink}")
    elif not all(isclose(incoming.get(name, 0), outgoing.get(name, 0), rel_tol=1e-6, abs_tol=1e-8)
                 for name in incoming.keys() | outgoing.keys()):
        drop = None
        diagnostics.append(f"Material utility {source}/{sink} has different component flows at its endpoints")
    return dict(inlet=source, outlet=sink, enthalpy_drop_kw=drop,
                heat_input_kw=None if drop is None else max(drop, 0.0))


def fuel_performance(model, plant, fuel_feeds, utility_pairs=()):
    """Product LHV divided by feed LHV and gross imported energy for this selection.

    Utility pairs are (inlet, outlet) IDs for the same circulating material.
    The fuel convention excludes sensible feed/product enthalpy, recovered
    heat and exported work credits. Missing constituent data propagate to the
    affected stream and efficiency while preserving the remaining results.
    """
    streams = {s["id"]: s for s in model["streams"]}
    components = {c["id"]: c for c in model["components"]}
    feeds = _selected(fuel_feeds, plant["inlets"], "Fuel feeds")
    products = _selected(plant["product_streams"], plant["outlets"], "Fuel products")
    if not products:
        raise ValueError("Fuel-performance accounting requires at least one selected product")
    if not isinstance(utility_pairs, (list, tuple)):
        raise ValueError("utility_pairs must contain (inlet, outlet) pairs")
    used = set(feeds + products)
    pairs, diagnostics = [], []
    for pair in utility_pairs:
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            raise ValueError("Expected a two-member (inlet, outlet) material utility pair")
        i, o = pair
        if i not in plant["inlets"] or o not in plant["outlets"] or i in used or o in used:
            raise ValueError(f"Material utility {i}/{o} overlaps another selection or uses the wrong boundary")
        pairs.append(_utility_pair(streams[i], streams[o], diagnostics))
        used.update(pair)

    species, fuel_streams = {}, {}
    for sid in feeds + products:
        rates = _rates(streams[sid])
        lhv_kw = None
        if rates is None:
            diagnostics.append(f"Missing or inconsistent component flows for {sid}")
        else:
            contributions = []
            for cid, n_kmol_h in rates.items():
                if n_kmol_h == 0:
                    continue
                if cid not in species:
                    if cid not in components:
                        raise ValueError(f"Unknown component {cid} on material stream {sid}")
                    species[cid] = _lhv_record(components[cid])
                lhv = species[cid]["lhv_j_mol"]
                if lhv is None:
                    diagnostics.append(f"Missing LHV for {cid} in {sid}: {species[cid]['reason']}")
                # kmol/h * J/mol converts to kW by division by 3600.
                contributions.append(None if lhv is None else n_kmol_h * lhv / 3600)
            lhv_kw = _total(contributions)
        fuel_streams[sid] = dict(role="feed" if sid in feeds else "product", lhv_kw=lhv_kw)
    feed_lhv, product_lhv = (_total(fuel_streams[s]["lhv_kw"] for s in side) for side in (feeds, products))
    material_heat = _total(pair["heat_input_kw"] for pair in pairs)
    heat_in, work_in = plant["heat_input_kw"], plant["work_input_kw"]
    denominator = _total((feed_lhv, heat_in, work_in, material_heat))
    efficiency = None
    if denominator is None:
        diagnostics.append("Gross energy input is incomplete")
    elif denominator <= 0:
        diagnostics.append("Fuel efficiency requires positive gross energy input")
    elif product_lhv is not None:
        efficiency = product_lhv / denominator
    return dict(
        basis="Selected product LHV / (selected feed LHV + gross heat in + gross work in + "
              "positive enthalpy drops of paired material utilities). No exported work or heat credits.",
        streams=fuel_streams, components=species, utility_pairs=pairs, diagnostics=diagnostics,
        product_lhv_kw=product_lhv, feed_lhv_kw=feed_lhv, purchased_heat_kw=heat_in,
        gross_work_input_kw=work_in, material_utility_heat_input_kw=material_heat,
        energy_input_kw=denominator, efficiency_fraction=efficiency,
    )
