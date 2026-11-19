"""Steady-flow control-volume balances. All transfer rates are positive inward."""

from collections import defaultdict
from math import fsum, isclose, isfinite, log, log1p


def _number(value, label):
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not isfinite(value):
        raise ValueError(f"{label} must be a finite number or null")
    return value


def _total(values):
    values = tuple(values)
    return None if None in values else fsum(values)


def _difference(incoming, outgoing):
    if incoming is None or outgoing is None:
        return None
    return incoming - outgoing


def _signed_parts(values):
    """Return the positive and negative magnitudes without losing missingness."""
    values = list(values)
    if None in values:
        return None, None
    return fsum(max(0, v) for v in values), fsum(max(0, -v) for v in values)


def _index(rows, name):
    index = {}
    for row in rows:
        key = row["id"]
        if not isinstance(key, str) or not key.strip() or key in index:
            raise ValueError(f"Invalid or duplicate {name} ID: {key!r}")
        index[key] = row
    return index


def _material_ports(blocks, streams):
    owners = {"inlets": {}, "outlets": {}}
    for block_id, block in blocks.items():
        for side, ports in owners.items():
            if not isinstance(block[side], (list, tuple)):
                raise ValueError(f"{block_id} {side} must be a list of material IDs")
            for stream_id in block[side]:
                if stream_id not in streams:
                    raise ValueError(f"{block_id} references unknown stream {stream_id}")
                if stream_id in ports:
                    raise ValueError(f"Multiple {side} use material stream {stream_id}")
                ports[stream_id] = block_id
    for stream_id, stream in streams.items():
        for key, side in (("from_block", "outlets"), ("to_block", "inlets")):
            owner = owners[side].get(stream_id)
            if key in stream and owner != stream[key]:
                raise ValueError(f"Stream {stream_id} {key} conflicts with block ports")
        if all(stream_id not in ports for ports in owners.values()):
            raise ValueError(f"Stream {stream_id} has no connection to the flowsheet")
    return owners["inlets"], owners["outlets"]


def _log_mean_temperature(low, high):
    """Effective boundary T for a constant heat-capacity rate."""
    if low == high:
        return low
    span = high - low
    ratio = span / low
    return span / (log1p(ratio) if isfinite(ratio) else log(high) - log(low))


def _transfer(block_id, transfer, t0_k):
    name = f"{block_id}/{transfer['id']}"
    kind = transfer["kind"]
    rate = _number(transfer["rate_kw"], f"{name} rate_kw")
    temperature = _number(transfer.get("temperature_k"), f"{name} temperature_k")
    limits = transfer.get("temperature_range_k")
    if limits is not None:
        if kind != "heat" or temperature is not None:
            raise ValueError(f"{name} temperature_range_k requires heat with no scalar temperature_k")
        if not isinstance(limits, (tuple, list)) or len(limits) != 2:
            raise ValueError(f"{name} temperature_range_k must have two temperatures in kelvin")
        points = [_number(t, f"{name} temperature_range_k") for t in limits]
        if any(t is None or t <= 0 for t in points):
            raise ValueError(f"{name} endpoints must be positive in kelvin")
        low, high = sorted(points)
        temperature = _log_mean_temperature(low, high)
        transfer = {**transfer, "temperature_k": temperature, "temperature_model": "constant_heat_capacity_rate"}
    if temperature is not None and temperature <= 0:
        raise ValueError(f"{name} boundary temperature must be positive in kelvin")
    if kind == "heat":
        if rate == 0:
            entropy = exergy = 0.0
        elif rate is None or temperature is None:
            entropy = exergy = None
        else:
            entropy = rate / temperature
            exergy = rate * (1.0 - t0_k / temperature)
    elif kind == "work":
        entropy = 0.0
        exergy = rate
    else:
        raise ValueError(f"{name} unknown transfer kind {kind!r}")
    if limits is not None and rate is not None and low < t0_k < high:
        # Heat crossing T0 carries exergy in opposite directions along its path.
        cold_heat = rate * ((t0_k - low) / (high - low))
        cold_exergy = cold_heat * (1 - t0_k / _log_mean_temperature(low, t0_k))
        ex_in, ex_out = _signed_parts((cold_exergy, exergy - cold_exergy))
    else:
        ex_in, ex_out = _signed_parts((exergy,))
    return dict(transfer, block=block_id, exergy_kw=exergy, entropy_kw_k=entropy,
                exergy_input_kw=ex_in, exergy_output_kw=ex_out)


def _internal_links(transfers):
    links = defaultdict(list)
    for leg in transfers:
        if leg.get("link_id") is not None:
            links[leg["link_id"]].append(leg)
    for name, legs in links.items():
        if len(legs) != 2 or legs[0]["block"] == legs[1]["block"]:
            raise ValueError(f"Internal link {name} needs one leg on each of two blocks")
        first, second = legs
        if first["kind"] != second["kind"]:
            raise ValueError(f"Internal link {name} has inconsistent transfer kinds")
        q1, q2 = first["rate_kw"], second["rate_kw"]
        if q1 is not None and q2 is not None and not isclose(q1, -q2, rel_tol=1e-6, abs_tol=1e-4):
            raise ValueError(f"Rates on internal link {name} do not cancel: {q1}, {q2} kW")
    return links


def _element_check(inlets, outlets, stream_results):
    """Account for atoms across reacting and nonreacting control volumes alike."""
    records = [stream_results[s].get("element_flows_kmol_h") for s in (*inlets, *outlets)]
    if any(record is None for record in records):
        return None, "incomplete"
    elements = sorted({element for record in records for element in record})
    residuals = {}
    status = "closed"
    for element in elements:
        flows = []
        for side in (inlets, outlets):
            flows.append(fsum(stream_results[s]["element_flows_kmol_h"].get(element, 0.0) for s in side))
        residuals[element] = flows[0] - flows[1]
        if not isclose(*flows, rel_tol=1e-6, abs_tol=1e-6):
            status = "residual"
    return residuals, status


def _boundary(name, inlets, outlets, streams, stream_results, transfers, t0_k):
    issues = []
    def stream_sum(ids, key, table):
        values = [_number(table[s].get(key), f"{s} {key}") for s in ids]
        missing = [s for s, v in zip(ids, values) if v is None]
        if missing:
            issues.append(f"Missing {key} for material stream(s) {', '.join(missing)}")
        return _total(values)

    mass_in, mass_out = (stream_sum(side, "mass_kg_h", streams) for side in (inlets, outlets))
    h_in, h_out = (stream_sum(side, "enthalpy_kw", streams) for side in (inlets, outlets))
    s_in, s_out = (stream_sum(side, "entropy_kw_k", streams) for side in (inlets, outlets))
    b_in, b_out = (stream_sum(side, "total_exergy_kw", stream_results) for side in (inlets, outlets))
    for leg in transfers:
        tag = f"{leg['block']}/{leg['id']}"
        if leg["rate_kw"] is None:
            issues.append(f"Missing {leg['kind']} rate for {tag}")
        elif leg["exergy_kw"] is None:
            issues.append(f"Missing boundary temperature for {tag}")

    q_in, q_out = _signed_parts(t["rate_kw"] for t in transfers if t["kind"] == "heat")
    w_in, w_out = _signed_parts(t["rate_kw"] for t in transfers if t["kind"] == "work")
    heat_ex, heat_entropy = (_total(t[key] for t in transfers if t["kind"] == "heat")
                             for key in ("exergy_kw", "entropy_kw_k"))
    transfer_ex_in, transfer_ex_out = (_total(t[key] for t in transfers)
                                       for key in ("exergy_input_kw", "exergy_output_kw"))
    ex_in, ex_out = _total((b_in, transfer_ex_in)), _total((b_out, transfer_ex_out))
    destruction = _difference(ex_in, ex_out)
    energy_in = _total((h_in, q_in, w_in))
    energy_out = _total((h_out, q_out, w_out))
    energy_residual = _difference(energy_in, energy_out)
    mass_residual = _difference(mass_in, mass_out)
    s_gen = _difference(s_out, _total((s_in, heat_entropy)))
    lost_work = None if s_gen is None else t0_k * s_gen

    # Scale closure checks by the largest independent term on the boundary.
    tolerance = None
    energy_status = "incomplete"
    if energy_residual is not None:
        tolerance = max(1e-3, 1e-7 * max(map(abs, (h_in, h_out, q_in, q_out, w_in, w_out))))
        energy_status = "closed" if abs(energy_residual) <= tolerance else "residual"
        if energy_status == "residual":
            issues.append("Energy residual exceeds the reported balance tolerance")
    if mass_residual is not None and not isclose(mass_in, mass_out, rel_tol=1e-6, abs_tol=1e-4):
        issues.append("Material mass inflow and outflow differ")
    exergy_status = "incomplete" if destruction is None else "available"
    if destruction is not None and destruction < -1e-3:
        exergy_status = "negative"
        issues.append("Negative signed exergy destruction exceeds 1 W")
    if lost_work is not None and lost_work < -1e-3:
        issues.append("Negative entropy generation exceeds 1 W on a T0 times Sgen basis")
    element_residuals, element_status = _element_check(inlets, outlets, stream_results)
    if element_status == "residual":
        issues.append("Atom inflow and outflow differ beyond elemental balance tolerance")
    reference_in, reference_out = (_total(stream_results[s].get("element_reference_kw") for s in side)
                                   for side in (inlets, outlets))
    reference_residual = _difference(reference_in, reference_out)
    predicted_balance = _total((lost_work, energy_residual, reference_residual))
    identity_residual = _difference(destruction, predicted_balance)
    if identity_residual is not None and abs(identity_residual) > max(1e-4, 1e-9 * max(abs(ex_in), abs(ex_out))):
        issues.append("Exergy balance differs from T0 Sgen plus energy and elemental-reference residuals")
    return {
        "id": name, "inlets": list(inlets), "outlets": list(outlets),
        "mass_in_kg_h": mass_in, "mass_out_kg_h": mass_out, "mass_residual_kg_h": mass_residual,
        "element_residual_kmol_h": element_residuals, "element_balance_status": element_status,
        "element_reference_residual_kw": reference_residual,
        "enthalpy_in_kw": h_in, "enthalpy_out_kw": h_out,
        "heat_input_kw": q_in, "heat_rejected_kw": q_out,
        "work_input_kw": w_in, "work_output_kw": w_out, "net_work_input_kw": _difference(w_in, w_out),
        "energy_residual_kw": energy_residual, "energy_balance_tolerance_kw": tolerance,
        "energy_balance_status": energy_status,
        "entropy_in_kw_k": s_in, "entropy_out_kw_k": s_out,
        "heat_entropy_net_kw_k": heat_entropy, "entropy_generation_kw_k": s_gen,
        "t0_entropy_generation_kw": lost_work, "heat_exergy_net_kw": heat_ex,
        "material_exergy_in_kw": b_in, "material_exergy_out_kw": b_out,
        "exergy_input_kw": ex_in, "exergy_output_kw": ex_out,
        "exergy_destruction_kw": destruction, "exergy_balance_status": exergy_status,
        "exergy_minus_t0_sgen_kw": _difference(destruction, lost_work),
        "exergy_identity_residual_kw": identity_residual,
        "diagnostics": issues,
    }


def _group_boundary(name, members, receivers, senders, transfers, links, results, streams, stream_results, t0_k):
    included = set(members)
    receiving = {s for s, block in receivers.items() if block in included}
    sending = {s for s, block in senders.items() if block in included}
    enclosed = {link for link, legs in links.items() if all(leg["block"] in included for leg in legs)}
    crossing = [leg for leg in transfers if leg["block"] in included and leg.get("link_id") not in enclosed]
    boundary = _boundary(name, sorted(receiving - sending), sorted(sending - receiving),
                         streams, stream_results, crossing, t0_k)
    boundary["members"] = list(members)
    boundary["internal_material_streams"] = sorted(sending & receiving)
    boundary["internal_energy_links"] = sorted(enclosed)
    boundary["crossing_transfers"] = crossing
    member_exergy = _total(results[bid]["exergy_destruction_kw"] for bid in members)
    member_energy = _total(results[bid]["energy_residual_kw"] for bid in members)
    # The thermal gradient on an enclosed link can account for additional irreversibility.
    link_exchange = _total(leg["exergy_kw"] for link in enclosed for leg in links[link])
    link_destruction = _difference(0.0, link_exchange)
    boundary["block_exergy_destruction_sum_kw"] = member_exergy
    boundary["internal_link_exergy_destruction_kw"] = link_destruction
    boundary["exergy_additivity_residual_kw"] = _difference(
        boundary["exergy_destruction_kw"], _total((member_exergy, link_destruction)))
    boundary["energy_additivity_residual_kw"] = _difference(boundary["energy_residual_kw"], member_energy)
    return boundary


def calculate_balances(model, stream_results, *, t0_k=298.15, product_streams=None, groups=None):
    """Balance blocks, named groups and the external plant boundary on a common datum.

    stream_results is keyed by material stream ID. Missing quantities are None.
    A link_id designates an internal transfer with two opposite signed legs.
    A product selection reports its share of all incoming boundary exergy.
    Each optional group maps a subsystem name to its participating block IDs.
    Overlapping groups describe separate boundaries and should not be summed.
    """
    if _number(t0_k, "t0_k") is None or t0_k <= 0:
        raise ValueError("t0_k must be positive in kelvin")
    blocks = _index(model["blocks"], "block")
    streams = _index(model["streams"], "stream")
    if not blocks:
        raise ValueError("The flowsheet must contain at least one block")
    if set(stream_results) != set(streams):
        raise ValueError("Stream exergy results and material inventory have different IDs")
    receivers, senders = _material_ports(blocks, streams)
    transfers = []
    results = {}
    for block_id, block in blocks.items():
        local = [_transfer(block_id, tr, t0_k) for tr in _index(block["transfers"], "transfer").values()]
        result = _boundary(block_id, block["inlets"], block["outlets"], streams, stream_results, local, t0_k)
        result["type"] = block.get("type")
        if "internal_heat_transfer_kw" in block:
            result["internal_heat_transfer_kw"] = _number(block["internal_heat_transfer_kw"], "internal_heat_transfer_kw")
        results[block_id] = result
        transfers.extend(local)

    links = _internal_links(transfers)
    groups = {} if groups is None else groups
    if not isinstance(groups, dict):
        raise ValueError("groups must map subsystem names to lists of participating block IDs")
    for name, members in groups.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("A subsystem group needs a nonempty string name")
        if (not isinstance(members, (list, tuple)) or not members
                or any(not isinstance(bid, str) or bid not in blocks for bid in members)):
            raise ValueError(f"Group {name} requires a nonempty list of known block IDs")
        if len(members) != len(set(members)):
            raise ValueError(f"Group {name} contains a duplicate member")

    plant, *subsystems = [
        _group_boundary(name, members, receivers, senders, transfers, links, results, streams, stream_results, t0_k)
        for name, members in [("plant", list(blocks)), *groups.items()]
    ]
    if isinstance(product_streams, str):
        raise ValueError("product_streams must be a sequence of external outlet IDs")
    products = list(product_streams) if product_streams is not None else []
    if len(set(products)) != len(products) or not set(products).issubset(plant["outlets"]):
        raise ValueError("Each product must be a distinct external material outlet")
    product_exergy = _total(stream_results[s].get("total_exergy_kw") for s in products) if products else None
    plant["product_streams"] = products
    plant["product_exergy_kw"] = product_exergy
    incoming = plant["exergy_input_kw"]
    fraction = None
    if product_exergy is not None and incoming is not None and incoming > 0:
        fraction = product_exergy / incoming
    plant["product_exergy_fraction_of_input"] = fraction
    return {"blocks": list(results.values()), "plant": plant, "transfers": transfers,
            "groups": {group["id"]: group for group in subsystems}}
