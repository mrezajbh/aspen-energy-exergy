"""Material energy and exergy at 298.15 K and 101325 Pa (rates in kW).

The elemental reference follows Tang et al., doi:10.3303/CET2183089, eq. 5.
The Aspen total requires enthalpy AND entropy referred to elemental formation.
"""

from collections import defaultdict
from math import fsum, isclose, isfinite, log

from ._estimate import PropertyUnavailable, estimate_physical, is_graphite
from ._exchange import case_hashes, digest

T0_K = 298.15
P0_PA = 101325.0
# J/mol of atoms. Tang et al., Table 1. Diatomic standards divided by two.
ELEMENT_EXERGY_J_MOL = {
    "C": 410260., "H": 118050., "N": 360., "O": 1985., "S": 609600.,
    "F": 233150., "Cl": 61800., "Br": 50600., "I": 9850., "Si": 854678.,
}
_REPORTED_QUALITIES = {"matched", "native", "hybrid", "reference"}
REFERENCE_PLAN_SCHEMA = "aspen-ambient-flash-plan/v1"
_REFERENCE_FIELDS = {
    "specific_SI": ("h0_J_kg", "s0_J_kgK"),
    "molar_SI": ("h0_J_mol", "s0_J_molK"),
    "flow_SI": ("H0_kW", "S0_kW_K"),
}


def _finite(value, label):
    if type(value) not in {float, int} or not isfinite(value):
        raise ValueError(f"{label} must be a finite number.")
    return value


def _nonnegative(value, label):
    if _finite(value, label) < 0:
        raise ValueError(f"{label} cannot be negative.")
    return value


def _fraction_vector(rates, total):
    if not isinstance(rates, dict) or type(total) not in (int, float) or not isfinite(total) or total <= 0:
        return None
    if any(type(n) not in (int, float) or not isfinite(n) or n < 0 for n in rates.values()):
        return None
    return {cid: n / total for cid, n in rates.items() if n > 0}


def _same_fractions(first, second):
    return (first is not None and second is not None and first.keys() == second.keys()
            and all(isclose(x, second[cid], rel_tol=1e-11, abs_tol=0.) for cid, x in first.items()))


def _reference_profile(stream):
    """Return fractional flows keyed by substream and component, on the whole-stream basis."""
    n = stream.get("molar_kmol_h")
    bulk = _fraction_vector(stream.get("component_flows_kmol_h"), n)
    subs = stream.get("substreams")
    if not bulk or not isclose(fsum(bulk.values()), 1., rel_tol=1e-8) or not isinstance(subs, dict) or not subs:
        return None
    profile = {}
    recovered = defaultdict(list)
    for sub, data in subs.items():
        if not isinstance(data, dict):
            return None
        fractions = _fraction_vector(data.get("component_flows_kmol_h"), n)
        if fractions is None:
            return None
        subflow = data.get("molar_kmol_h")
        if subflow is not None and (type(subflow) not in (int, float) or not isfinite(subflow)
                                    or not isclose(fsum(fractions.values()), subflow / n, rel_tol=1e-8)):
            return None
        for cid, frac in fractions.items():
            profile[sub, cid] = frac
            recovered[cid].append(frac)
    if not _same_fractions(bulk, {cid: fsum(xs) for cid, xs in recovered.items()}):
        return None
    return profile


def equivalent_reference_streams(first, second):
    """Check saved evidence for matching reference-flash feeds, independent of throughput.

    Operating phases may differ across temperatures. Ambient phase permissions
    and hidden property options additionally require the exporter's confirmation.
    """
    return (bool(first.get("property_method")) and first.get("property_method") == second.get("property_method")
            and first.get("property_options_id") == second.get("property_options_id")
            and first.get("reference_options_key") == second.get("reference_options_key")
            and _same_fractions(_reference_profile(first), _reference_profile(second)))


def _check_plan_feed(record, stream):
    sid = stream["id"]
    if "operating_state_sha256" in record and record["operating_state_sha256"] != digest(stream):
        raise ValueError(f"{sid}: operating-state fingerprint differs from the reference plan.")
    if "reference_options_key" in record and record["reference_options_key"] != stream.get("reference_options_key"):
        raise ValueError(f"{sid}: saved property options differ from the reference plan.")
    if "substream_component_flows_kmol_h" not in record:
        return
    feeds = record["substream_component_flows_kmol_h"]
    subs = stream.get("substreams") or {}
    if not isinstance(feeds, dict) or feeds.keys() != subs.keys():
        raise ValueError(f"{sid}: reference plan substream names differ from the case.")
    expected = _reference_profile(stream)
    if expected is None:
        if feeds != {k: v.get("component_flows_kmol_h") for k, v in subs.items()}:
            raise ValueError(f"{sid}: incomplete substream inventory differs from the case.")
        return
    flows = {}
    for sub, rates in feeds.items():
        if not isinstance(rates, dict):
            raise ValueError(f"{sid}: missing component flows in reference plan substream {sub}.")
        for cid, rate in rates.items():
            flows[sub, cid] = _nonnegative(rate, f"{sid}/{sub}/{cid} reference feed")
    if not _same_fractions(expected, _fraction_vector(flows, fsum(flows.values()))):
        raise ValueError(f"{sid}: reference plan substream component allocation differs from the case.")


def _reference_index(document, model):
    if document is None:
        return {}
    if isinstance(document, (list, tuple)):
        combined = {}
        for part in document:
            accepted = _reference_index(part, model)
            duplicate = combined.keys() & accepted.keys()
            if duplicate:
                raise ValueError(f"Multiple reference returns cover the same stream(s): {', '.join(sorted(duplicate))}")
            combined.update(accepted)
        return combined
    if not isinstance(document, dict):
        raise ValueError("Ambient reference document must be a JSON object.")
    case_digest = document.get("case_sha256")
    if not isinstance(case_digest, str) or case_digest.lower() not in case_hashes(model):
        raise ValueError("Ambient reference document SHA-256 does not match this Aspen case.")
    marker = document.get("plan_schema")
    if marker is not None and marker != REFERENCE_PLAN_SCHEMA:
        raise ValueError(f"Unknown ambient reference plan schema {marker!r}.")
    records = document.get("streams")
    if not isinstance(records, list):
        raise ValueError("Ambient reference document requires a streams list.")
    streams = {s["id"]: s for s in model["streams"]}
    indexed, representatives = {}, set()
    for record in records:
        sid = record.get("id") if isinstance(record, dict) else None
        if not isinstance(sid, str) or sid not in streams or sid in representatives:
            raise ValueError(f"Duplicate or unknown ambient reference stream {sid!r}.")
        representatives.add(sid)
        if marker == REFERENCE_PLAN_SCHEMA:
            return_keys = ["aspen_reference_id", "source_location",
                           *(field for pair in _REFERENCE_FIELDS.values() for field in pair)]
            if (all(record.get(k) in (None, "") for k in return_keys)
                    and record.get("same_aspen_property_options") is not True
                    and record.get("same_aspen_property_options_for_all") is not True):
                continue
            _check_plan_feed(record, streams[sid])
        peers = record.get("applies_to", [])
        if not isinstance(peers, list) or any(not isinstance(peer, str) for peer in peers):
            raise ValueError(f"{sid}: applies_to must list other operating stream IDs.")
        if peers and record.get("same_aspen_property_options_for_all") is not True:
            raise ValueError(f"{sid}: confirm common Aspen property, equilibrium-phase and substream options "
                             "for all applies_to streams.")
        for target in [sid, *peers]:
            if target not in streams or target in indexed:
                raise ValueError(f"Duplicate or unknown ambient reference target {target!r}.")
            if target != sid:
                origin, peer = streams[sid], streams[target]
                if not equivalent_reference_streams(origin, peer):
                    raise ValueError(f"{sid} and {target}: cannot verify common composition, substreams or property options.")
                if record.get("quantity_basis") == "specific_SI":
                    mw = _nonnegative(origin.get("mass_kg_h"), f"{sid} mass") / origin["molar_kmol_h"]
                    peer_mw = _nonnegative(peer.get("mass_kg_h"), f"{target} mass") / peer["molar_kmol_h"]
                    if mw <= 0 or not isclose(mw, peer_mw, rel_tol=1e-8, abs_tol=0.):
                        raise ValueError(f"{sid} and {target}: mass-specific reference requires matching molar mass.")
            indexed[target] = record
    return indexed


def _reference_availability(record, stream, component_ids):
    """Validate source, method, equilibrium, reference conditions and composition."""
    sid = stream["id"]
    for key in ("aspen_reference_id", "source_location"):
        if not isinstance(record.get(key), str) or not record[key].strip():
            raise ValueError(f"{sid}: missing reference provenance field {key}.")
    if (record.get("phase_convention") != "equilibrium"
            or record.get("same_aspen_property_options") is not True
            or not stream.get("property_method")
            or record.get("property_method") != stream["property_method"]):
        raise ValueError(f"{sid}: reference requires an equilibrium flash with the same Aspen property method and options.")
    for key, target, tol in (("T_K", T0_K, 1e-6), ("P_Pa", P0_PA, 1e-4)):
        if not isclose(_finite(record.get(key), key), target, abs_tol=tol, rel_tol=0):
            raise ValueError(f"{sid}: incompatible ambient reference {key}.")
    flows = record.get("component_flows_kmol_h")
    if record.get("component_ids") != component_ids or not isinstance(flows, list) or len(flows) != len(component_ids):
        raise ValueError(f"{sid}: reference component list and order must match the Aspen inventory.")
    ref_n = fsum(_nonnegative(n, f"{sid} reference molar flow") for n in flows)
    n = stream["molar_kmol_h"]
    if not ref_n or not n:
        raise ValueError(f"{sid}: reference composition requires positive flow.")
    actual = stream["component_flows_kmol_h"]
    if record.get("applies_to"):
        composition_matches = _same_fractions(_fraction_vector(actual, n),
                                              _fraction_vector(dict(zip(component_ids, flows)), ref_n))
    else:
        composition_matches = all(isclose(actual.get(cid, 0.) / n, v / ref_n, rel_tol=1e-9, abs_tol=1e-12)
                                  for cid, v in zip(component_ids, flows))
    if not composition_matches:
        raise ValueError(f"{sid}: ambient and operating compositions differ.")
    fields = _REFERENCE_FIELDS
    basis = record.get("quantity_basis")
    if basis not in fields:
        raise ValueError(f"{sid}: unknown reference quantity_basis.")
    keys = fields[basis]
    if any(k in record for pair in fields.values() for k in pair if k not in keys):
        raise ValueError(f"{sid}: conflicting reference property bases.")
    h0, s0 = (_finite(record.get(k), k) for k in keys)
    if basis == "specific_SI":
        scale = _nonnegative(stream.get("mass_kg_h"), f"{sid} operating mass") / 3.6e6
    elif basis == "molar_SI":
        scale = n / 3600.
    else:
        scale = n / ref_n
    return (h0 - T0_K * s0) * scale


def _composition(stream, components):
    sid = stream["id"]
    n = stream.get("molar_kmol_h")
    if n is not None:
        _nonnegative(n, f"{sid} molar flow")
    mass = stream.get("mass_kg_h")
    if mass is not None:
        _nonnegative(mass, f"{sid} mass flow")
    if n == 0 and mass not in (None, 0):
        raise ValueError(f"{sid}: zero molar flow with nonzero mass flow.")
    if mass == 0 and n is not None and n > 0:
        raise ValueError(f"{sid}: zero mass flow with positive molar flow.")
    flows = stream.get("component_flows_kmol_h")
    if flows is None:
        return [] if n == 0 else None
    if not isinstance(flows, dict):
        raise ValueError(f"{sid}: component_flows_kmol_h must map component IDs to rates.")
    if unknown := flows.keys() - components.keys():
        raise ValueError(f"{sid}: unrecognized components {sorted(unknown)}.")
    for cid, rate in flows.items():
        _nonnegative(rate, f"{sid} component flow {cid}")
    if n == 0 and any(flows.values()):
        raise ValueError(f"{sid}: component flow in an empty stream.")
    if n is not None and not isclose(fsum(flows.values()), n, abs_tol=1e-9, rel_tol=1e-8):
        raise ValueError(f"{sid}: component rates do not sum to the total molar flow.")
    if not any(flows.values()) and n != 0:
        return None
    return [(components[cid], rate) for cid, rate in flows.items() if rate > 0]


def _elements(active):
    contributions = {}
    for component, rate in active:
        atoms = component.get("atoms")
        if not atoms:
            return None, None
        for element, count in atoms.items():
            _nonnegative(count, f"{component['id']} atom count {element}")
            contributions.setdefault(element, []).append(rate * count)
    flows = {el: fsum(terms) for el, terms in contributions.items()}
    if flows.keys() - ELEMENT_EXERGY_J_MOL.keys():
        return flows, None
    term = fsum(rate * ELEMENT_EXERGY_J_MOL[el] for el, rate in flows.items()) / 3600.
    return flows, term


def _stream(stream, components, formation, reference, estimate_splits, native):
    sid = stream["id"]
    h, s = stream.get("enthalpy_kw"), stream.get("entropy_kw_k")
    for key in ("enthalpy_kw", "entropy_kw_k", "temperature_k", "pressure_pa", "vapor_fraction"):
        val = stream.get(key)
        if val is not None:
            _finite(val, f"{sid} {key}")
    vf = stream.get("vapor_fraction")
    if vf is not None and not -1e-10 <= vf <= 1 + 1e-10:
        raise ValueError(f"{sid}: vapor fraction outside [0, 1].")
    active = _composition(stream, components)
    atomic_flows, term = _elements(active) if active is not None else (None, None)
    notes = []
    total = physical = chemical = None
    availability = None if h is None or s is None else h - T0_K * s
    method, quality = "unavailable", "unavailable"
    if active == []:
        if h not in (None, 0) or s not in (None, 0):
            raise ValueError(f"{sid}: empty stream carries nonzero enthalpy or entropy.")
        if reference is not None:
            raise ValueError(f"{sid}: ambient reference supplied for a zero-flow stream.")
        total = physical = chemical = 0.
        method, quality = "zero_flow", "zero"
    else:
        if formation and term is not None and availability is not None:
            total = fsum((h, -T0_K * s, term))
        else:
            if not formation:
                notes.append("Elemental-formation enthalpy/entropy datum has not been confirmed.")
            if term is None:
                notes.append("Missing elemental composition or environmental standard for a flowing component.")
            if availability is None:
                notes.append("Missing operating enthalpy or entropy rate.")
        if reference is not None:
            if active is None or stream.get("molar_kmol_h") is None:
                raise ValueError(f"{sid}: ambient reference requires operating composition and molar flow.")
            ref = _reference_availability(reference, stream, list(components))
            if availability is not None:
                physical = fsum((h, -T0_K * s, -ref))
            if formation and term is not None:
                chemical = ref + term
            method, quality = "aspen_ambient_flash", "matched"
            notes.append(f"{reference['aspen_reference_id']} ({reference['source_location']})")
            if sid != reference["id"]:
                method = "aspen_shared_ambient_flash"
                notes.append(f"Reference flash for {reference['id']} scaled to {sid} after confirmation "
                             "of common Aspen property, phase and substream options.")
        elif total is not None and len(active) == 1 and is_graphite(active[0][0]):
            physical, chemical = availability, term
            method, quality = "elemental_graphite", "reference"
        elif native is None and active and estimate_splits:
            try:
                physical, method, model_notes = estimate_physical(stream, active, T0_K, P0_PA, formation=formation)
                chemical = None if total is None else total - physical
                quality = "estimated"
                notes.extend(model_notes)
                if {"Liquid1", "Liquid2"}.issubset(stream.get("active_phases") or []):
                    notes.append("Aspen reports two liquid phases. The independent model permits only one.")
                notes.append("Independent physical exergy uses paired model states. Its chemical counterpart is the Aspen total minus that estimate.")
            except PropertyUnavailable as exc:
                notes.append(str(exc))
        elif native is None:
            notes.append("Same-composition ambient Aspen flash required for this split.")
    if native is not None:
        reported = _finite(native.get("physical_exergy_kw"), f"{sid} native EXERGYFL")
        if active == [] and reported != 0:
            raise ValueError(f"{sid}: nonzero Aspen EXERGYFL for an empty stream.")
        if physical is not None and not isclose(reported, physical, rel_tol=1e-6, abs_tol=1e-3):
            raise ValueError(f"{sid}: Aspen EXERGYFL ({reported:.9g} kW) disagrees with "
                             f"{method} ({physical:.9g} kW). Check the reference and substream scope.")
        physical = reported
        if total is not None:
            chemical = total - physical
        hybrid = native["scope"] == "mixed_with_saved_graphite"
        method = "aspen_exergyfl_with_saved_graphite" if hybrid else "aspen_exergyfl"
        quality = "hybrid" if hybrid else "native"
        notes.append(f"Aspen EXERGYFL source: {native['export_source']}. Scope: {native['scope']}.")
        if hybrid:
            notes.append("Physical exergy includes the formation-referenced saved graphite contribution.")
    provenance = {} if reference is None else {
        "reference_stream_id": reference["id"], "aspen_reference_id": reference["aspen_reference_id"],
        "reference_source_location": reference["source_location"], "reference_quantity_basis": reference["quantity_basis"],
    }
    return {
        "id": sid, "enthalpy_kw": h, "entropy_kw_k": s, "element_flows_kmol_h": atomic_flows,
        "element_reference_kw": term, "total_exergy_kw": total, "physical_exergy_kw": physical,
        "chemical_exergy_kw": chemical, "split_method": method, "split_quality": quality, "notes": notes,
        **provenance, **({"native_exergy": dict(native)} if native is not None else {}),
    }


def _share_ambient_estimates(streams, rows, *, estimate_splits):
    """Use one ambient chemical potential for equivalent source compositions.

    Save the original two-state estimate for comparison. Native enthalpy and
    entropy differences then carry changes between two equivalent streams.
    """
    groups = defaultdict(list)
    for stream in streams:
        n = stream.get("molar_kmol_h")
        row = rows[stream["id"]]
        if not n or row["total_exergy_kw"] is None or not stream.get("property_method"):
            continue
        if stream.get("substreams"):
            fractions = _reference_profile(stream)
        else:
            # Normalized exports can omit substream data. In this mode shared
            # properties remain explicitly estimated unless a flash is returned.
            fractions = _fraction_vector(stream["component_flows_kmol_h"], n)
        if not fractions:
            continue
        signature = (stream["property_method"], tuple(sorted(stream.get("substreams") or {})),
                     stream.get("property_options_id"), stream.get("reference_options_key"), tuple(sorted(fractions)))
        candidates = groups[signature]
        group = next((members for known, members in candidates if _same_fractions(fractions, known)), None)
        if group is None:
            group = []
            candidates.append((fractions, group))
        group.append(stream)
    for candidates in groups.values():
        for _, members in candidates:
            if len(members) < 2:
                continue
            def rank(stream):
                row = rows[stream["id"]]
                t, p = stream.get("temperature_k"), stream.get("pressure_pa")
                distance = abs(log(t / T0_K)) + abs(log(p / P0_PA)) if t and p and min(t, p) > 0 else float("inf")
                # At nearly ambient T and P, prefer the smallest modeled molar
                # correction. A few kPa of liquid pressure can cost much less
                # availability than a few degrees of temperature change.
                near_ambient = t and p and abs(t / T0_K - 1) <= 0.05 and abs(p / P0_PA - 1) <= 0.15
                correction = abs(row["physical_exergy_kw"] or 0.) / stream["molar_kmol_h"]
                return (row["split_quality"] not in _REPORTED_QUALITIES, row["split_quality"] == "hybrid", not near_ambient,
                        correction if near_ambient else distance, distance, stream["id"])
            choices = [s for s in members if rows[s["id"]]["chemical_exergy_kw"] is not None]
            if not choices:
                continue
            anchor = min(choices, key=rank)
            source = rows[anchor["id"]]
            chemical_per_kmol_h = source["chemical_exergy_kw"] / anchor["molar_kmol_h"]
            conflicts = []
            for stream in choices:
                row = rows[stream["id"]]
                if row["split_quality"] not in _REPORTED_QUALITIES:
                    continue
                expected = chemical_per_kmol_h * stream["molar_kmol_h"]
                if not isclose(row["chemical_exergy_kw"], expected, rel_tol=1e-6, abs_tol=1e-3):
                    conflicts.append(stream["id"])
                    row["comparison_reference_id"] = anchor["id"]
                    row["reference_difference_kw"] = row["chemical_exergy_kw"] - expected
            if conflicts:
                ids = ", ".join([anchor["id"], *conflicts])
                for peer in members:
                    rows[peer["id"]]["notes"].append(
                        f"Reference streams {ids} imply different molar chemical exergies. "
                        "This composition group retains each stream's individual result.")
                continue
            if not estimate_splits:
                continue
            for peer in members:
                row = rows[peer["id"]]
                if row is source or row["split_quality"] in _REPORTED_QUALITIES:
                    continue
                old_physical = row["physical_exergy_kw"]
                old_method = row["split_method"]
                chemical = source["chemical_exergy_kw"] * (peer["molar_kmol_h"] / anchor["molar_kmol_h"])
                physical = row["total_exergy_kw"] - chemical
                row.update(chemical_exergy_kw=chemical, physical_exergy_kw=physical, split_quality="estimated",
                           split_method="shared_ambient_reference", reference_anchor_id=anchor["id"],
                           reference_anchor_method=source["split_method"], independent_split_method=old_method,
                           independent_physical_exergy_kw=old_physical,
                           independent_difference_kw=None if old_physical is None else physical - old_physical)
                row["notes"].append(
                    f"Reported split uses the molar chemical exergy of {anchor['id']} at matching composition and "
                    "Aspen method, assuming common equilibrium ambient options. Native availability differences "
                    "give the reported physical exergy. The independent property estimate is retained separately.")


def calculate_streams(model, *, references=None, estimate_splits=True, native_exergy=None):
    """Compute total flow exergy and estimate or match the physical/chemical split.

    Pass the optional ambient-reference document as a dict and validated Aspen
    EXERGYFL records from read_native_exergy as native_exergy. Disable model
    estimates to keep supplied Aspen values and the elemental graphite standard.
    Chemical standards apply to the listed elements at T0_K and P0_PA.
    """
    components = {c["id"]: c for c in model["components"]}
    metadata = model.get("metadata", {})
    formation = metadata.get("enthalpy_entropy_basis", metadata.get("property_datum")) == "aspen_formation"
    streams = model["streams"]
    stream_ids = {s["id"] for s in streams}
    refs = _reference_index(references, model)
    native = {} if native_exergy is None else native_exergy
    if not isinstance(native, dict) or native.keys() - stream_ids:
        raise ValueError("Native EXERGYFL records must be indexed by known stream IDs.")
    rows = {s["id"]: _stream(s, components, formation, refs.get(s["id"]), estimate_splits,
                              native.get(s["id"])) for s in streams}
    if estimate_splits or native or refs:
        _share_ambient_estimates(streams, rows, estimate_splits=estimate_splits)
    return rows

