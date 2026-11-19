"""Read solved Aspen Plus 41.0 BKP/APWZ caches without starting Aspen.

Properties are taken from the embedded ADS XML in SI. The text flowsheet and
stream result headers supply the port connections, phase labels and methods.
"""

import base64
from collections import defaultdict
import hashlib
from io import BytesIO
import json
from math import fsum, isclose, isfinite
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile
import zlib


# Local integer identifiers can change when Aspen serializes another case.
# Resolve by property UUID instead of assuming positions or component counts.
_UUID = {
    "block_id": "4e87a2bd-421b-4b63-8ed8-79e387678252",
    "heat_w": "2348d9d1-decf-4beb-b6cf-571f5f38f08c",
    "material_record": "ebe609c5-c41b-4f79-88a6-ed3c0abfb958",
    "stream_id": "ba430646-cca5-47e0-b6f0-7051e3417088",
    "substreams": "c690f9dd-0a30-4ae6-8055-34c0fc58c300",
    "temperature": "36adff82-c12d-4b22-9ac5-1e4bc786ed91",
    "pressure": "aac6dd64-f01f-4ebf-b4c8-57c4b850a9be",
    "mass_flow": "4ed0cdbb-bb83-4e2d-a59d-5adf243aef69",
    "molar_flow": "d4198184-cf6e-41e8-b963-e3b5b1352ebb",
    "vapor_fraction": "5729a872-3fa6-47ad-b309-341a16bcaa62",
    "liquid_fraction": "f47d9357-2c29-4fff-ad60-e9a0e39a2fe2",
    "solid_fraction": "c3ecab21-3cac-45c0-9709-8ad2313b8457",
    "h_j_kg": "643f1978-21c9-4529-a95c-3f8e52072309",
    "s_j_kg_k": "3fc900c7-d0f0-4f28-9799-a033461dbab2",
    "mole_fractions": "8e4f9fb5-1f4a-4a9c-bc51-e561eddf40b7",
    "vapor_record": "2194afbf-4c01-424a-a598-2d6b41502805",
    "liquid1_record": "565a9ba8-597c-4394-9361-89564924fe3e",
    "liquid2_record": "9f0bd91f-1b59-4473-9f30-cccf63f6787f",
    "phase_mass_flow": "46a17462-bfda-4b06-8d89-caa76c57cf7c",
    "component_record": "d2a1d229-925a-4ce5-bf56-f766fc3e263b",
    "component_id": "d02bfd0b-3110-4c72-9316-ae824f8c0a4f",
    "alias": "2431fc50-91f7-440a-a87c-f0502f54b68e",
    "mw": "20a61041-f1e7-42dd-9079-7e90a1d64c32",
    "compressor_shaft_w": "616c0999-5ffc-4968-ad03-a0c4919c9284",
    "compressor_loss_w": "1d377e6a-992f-43a3-9886-028042eb9061",
    "pump_fluid_w": "fe83efee-d632-4f03-873f-ddc318a7dfd3",
    "pump_shaft_w": "ea50fef7-d542-4d0c-a0c8-e31f3e980033",
    "condenser_w": "1ef4a839-7b5b-4f38-b654-3657cd39c78a",
    "reboiler_w": "b68d877e-e5db-4d25-baaf-dbca562adc93",
    "condenser_k": "60ab53c9-5723-4fa9-84df-bb7e8f6a06b0",
    "reboiler_k": "73a14b3d-546e-42e4-8ff2-cc33f2d83c80",
}
_AUDITED_FORMATION_BACKUP = "4e16eb57942f8db308ee6ac44863b7301eeaa60704b3ce2cea2809a40fc45fc2"
EXERGY_REPORT_BACKUP = "b188fa23b2e33f58f8a85db9159a4faa8d0f70d56bc0a794f22d2c3b41719375"
_ID = r'("(?:[^"]|"")*"|[^\s()]+)'
_SCALAR = re.compile(r'"(?:[^"]|"")*"|[^\s()]+')
_HEAT_MODELS = {"HEATER", "SEP", "FLASH2", "RSTOIC", "RGIBBS", "RYIELD"}
_PASSIVE_MODELS = {"MIXER", "FSPLIT", "SSPLIT", "VALVE"}


def _name(token):
    return token[1:-1].replace('""', '"') if token.startswith('"') else token


def _unique(mapping, key, value):
    if key in mapping:
        raise ValueError(f"Duplicate Aspen record {key!r}.")
    mapping[key] = value


def _number(value, field):
    try:
        result = float(value.replace("D", "E"))
    except (AttributeError, ValueError) as error:
        raise ValueError(f"Missing or invalid saved Aspen field {field}.") from error
    if not isfinite(result):
        raise ValueError(f"Non-finite saved Aspen field {field}.")
    return result


def _backup(path, source):
    if path.suffix.lower() == ".bkp":
        return path.name, source
    if path.suffix.lower() != ".apwz":
        raise ValueError("Expected a saved Aspen Plus .bkp or .apwz file.")
    try:
        with zipfile.ZipFile(BytesIO(source)) as archive:
            members = [n for n in archive.namelist() if n.lower().endswith(".bkp")]
            if len(members) != 1:
                raise ValueError(f"Expected one .bkp inside the APWZ, found {len(members)}.")
            return members[0], archive.read(members[0])
    except zipfile.BadZipFile as error:
        raise ValueError("Invalid Aspen APWZ ZIP archive.") from error


def _ports(group):
    tokens = [_name(s) for s in _SCALAR.findall(group)]
    if len(tokens) % 2:
        raise ValueError("Incomplete Aspen flowsheet port pair.")
    material, energy = [], []
    for sid, port in zip(tokens[::2], tokens[1::2]):
        if port.startswith("M"):
            material.append(sid)
        elif port.startswith(("Q", "W")):
            energy.append(sid)
        else:
            raise ValueError(f"Unsupported Aspen port {port} for {sid}.")
    return material, energy


def _topology(flat):
    match = re.search(r"FLOWSHEET GLOBAL\b(.*?)\bPROPERTIES MAIN\b", flat)
    if not match:
        raise ValueError("Missing supported Aspen FLOWSHEET GLOBAL section.")
    flow = match[1]
    block_pattern = (r"BLOCK BLKID\s*=\s*" + _ID + r" BLKTYPE\s*=\s*" + _ID
                     + r" MDLTYPE\s*=\s*" + _ID + r" IN\s*=\s*\((.*?)\) OUT\s*=\s*\((.*?)\)")
    blocks = {}
    for row in re.finditer(block_pattern, flow):
        bid, kind = _name(row[1]), _name(row[2])
        ins, ein = _ports(row[4])
        outs, eout = _ports(row[5])
        _unique(blocks, bid, dict(id=bid, type=kind, inlets=ins, outlets=outs,
                                 transfers=[], energy_inlets=ein, energy_outlets=eout))
    if not blocks or len(blocks) != len(re.findall(r"\bBLOCK BLKID\s*=", flow)):
        raise ValueError("Incomplete or unsupported Aspen block port layout.")
    return blocks


def _result_headers(text):
    streams, links, warnings = {}, {}, []
    header = r"^ DSET STREAM (\w+) " + _ID + r" STR_MAIN @(\S+)\s*\("
    for m in re.finditer(header, text, re.MULTILINE):
        kind, sid, layout = m[1], _name(m[2]), m[3]
        # Read only the header. ADS supplies all material property values.
        parts = _STREAM_FIELDS.match(text, m.end())
        if parts is None:
            raise ValueError(f"Incomplete {sid} STR_MAIN result header.")
        origin, target, phase, fourth = (_name(x) for x in parts.groups())
        end = dict(id=sid, from_block=None if origin == "*" else origin,
                   to_block=None if target == "*" else target)
        if kind == "MATERIAL":
            _unique(streams, sid, {**end, "phase": phase, "property_method": fourth})
        elif (kind, layout) in {("HEAT", "L_HEAT"), ("WORK", "L_WORK")}:
            unit = re.compile(r"\s*<(\d+)>").match(text, parts.end())
            expected = "15" if kind == "HEAT" else "3"
            if not unit or unit[1] != expected:
                raise ValueError(f"Unsupported {kind} display-unit code at {sid}.")
            rate = _number(fourth, sid) * (1000 if kind == "HEAT" else 1)
            _unique(links, sid, {**end, "kind": kind.lower(), "saved_rate_kw": rate})
        else:
            raise ValueError(f"Unsupported Aspen energy stream layout {sid} ({kind}/{layout}).")
    pattern = r"^ DSET (?:BLOCK|STREAM) \S+ " + _ID + r" ERRTXT @\S+\s*\((.*?)\)[ \t]*$"
    for m in re.finditer(pattern, text, re.MULTILINE | re.DOTALL):
        lines = [_name(t).strip() for t in re.findall(r'"(?:[^"]|"")*"', m[2])]
        warnings.append(f"{_name(m[1])}: " + " ".join(lines))
    return streams, links, warnings


_STREAM_FIELDS = re.compile(r"\s*" + _ID + r"\s+" + _ID + r"\s+" + _ID + r"\s+" + _ID)


def _check_topology(blocks, stream_headers, links):
    for records, sides in ((stream_headers, ("inlets", "outlets")),
                           (links, ("energy_inlets", "energy_outlets"))):
        for side, endpoint in zip(sides, ("to_block", "from_block")):
            actual = {}
            for b in blocks.values():
                for sid in b[side]:
                    _unique(actual, sid, b["id"])
            expected = {sid: row[endpoint] for sid, row in records.items() if row[endpoint] is not None}
            if actual != expected:
                raise ValueError(f"Aspen {side} differ between flowsheet and saved streams.")


def _component_info(flat):
    section = re.search(r"COMPONENTS MAIN\s*\?(.*?)(?:\?\s*COMPONENTS\s|$)", flat)
    info = {}
    if section:
        pattern = (r"\bCID\s*=\s*" + _ID + r" ANAME\s*=\s*" + _ID + r" OUTNAME\s*=\s*" + _ID
                   + r"(.*?)\bDBNAME1\s*=\s*" + _ID + r" ANAME1\s*=\s*" + _ID)
        for m in re.finditer(pattern, section[1]):
            formula = _name(m[6]).split("-")[0]
            # The Aspen database alias convention contains a formula followed
            # by an isomer suffix. Reject formulas requiring a different grammar.
            if not re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", formula):
                raise ValueError(f"Unsupported component formula {formula!r} for {_name(m[1])}.")
            atoms = defaultdict(int)
            for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
                atoms[element] += int(count or 1)
            _unique(info, _name(m[1]), dict(name=_name(m[5]), formula=formula, atoms=dict(atoms),
                                          type="solid" if "SOLID" in m[4] else "conventional"))
    return info


def _deck_parts(text, delimiter):
    """Split Aspen input records at an unquoted delimiter."""
    start = 0
    for token in re.finditer(r'"(?:[^"]|"")*"|' + re.escape(delimiter), text):
        if token[0] == delimiter:
            yield text[start:token.start()].strip()
            start = token.end()
    yield text[start:].strip()


def _reference_option_hints(flat, streams):
    """Record explicit property/phase settings without inferring Aspen defaults."""
    sections, section_options, block_records = {}, {}, {}
    deck = list(_deck_parts(flat, "?"))
    for header, body in zip(deck, deck[1:]):
        flow = re.fullmatch(r"FLOWSHEET " + _ID, header)
        block = re.fullmatch(r"BLOCK [\w-]+ " + _ID, header)
        if not (flow or block or header == "PROPERTIES MAIN"):
            continue
        commands = {}
        for command in _deck_parts(body, "\\"):
            tokens = _SCALAR.findall(command)
            if not tokens:
                continue
            name = _name(tokens[0])
            if name not in {"BLOCK", "GPROPERTIES", "PROPERTIES", "PARAM", "BLOCK-OPTION",
                            "OUTLETS-HOT", "OUTLETS-COLD"}:
                continue
            fields = {k: _name(v) for k, v in re.findall(r"(?:^|\s)([A-Z][A-Z0-9-]*)\s*=\s*" + _ID, command)}
            if flow and name == "BLOCK" and "BLKID" in fields:
                sections[fields["BLKID"]] = _name(flow[1])
            elif block:
                commands.setdefault(name, {}).update(fields)
            elif header == "PROPERTIES MAIN" and name == "GPROPERTIES":
                section_options["GLOBAL"] = fields
            elif header == "PROPERTIES MAIN" and name == "PROPERTIES" and "SECID-LIST" in fields:
                section_options[fields.pop("SECID-LIST")] = fields
        if block:
            block_records[_name(block[1])] = commands

    side_fields = ("OPSET", "FREE-WATER", "SOLU-WATER", "TRUE-COMPS")
    sided_keys = {prefix + key for prefix in ("C", "H") for key in side_fields}
    for stream in streams:
        producer = stream["from_block"]
        commands = block_records.get(producer, {})
        section = sections.get(producer)
        raw = commands.get("PROPERTIES", {})
        side = None
        for label, prefix in (("HOT", "H"), ("COLD", "C")):
            if commands.get("OUTLETS-" + label, {}).get(label + "-SID") == stream["id"]:
                side = prefix
        options = {k: v for k, v in raw.items() if k not in sided_keys}
        if side:
            options.update({k: raw[side + k] for k in side_fields if side + k in raw})
        phase = {k: v for k, v in commands.get("PARAM", {}).items() if k in {"NPHASE", "NO-PHASE", "PHASE"}}
        # GFREE-WATER/FREE-WATER select a water property method. Only the
        # explicit producer BLKOPFREWAT field establishes this block's switch.
        hints = dict(section=section, section_options=section_options.get(section, {}),
                     producer_free_water=commands.get("BLOCK-OPTION", {}).get("BLKOPFREWAT"),
                     producer_phase_options=phase, producer_property_options=options)
        if not side and raw.keys() & sided_keys:
            hints["unresolved_exchanger_outlet"] = stream["id"]
        key = json.dumps(hints, sort_keys=True, separators=(",", ":")).encode("utf-8")
        stream["reference_options_key"] = "aspen41:" + hashlib.sha256(key).hexdigest()
        stream["reference_option_hints"] = dict(hints, producer_block=producer, exchanger_side=side)
    return section_options


class _ADS:
    def __init__(self, tree):
        dictionary, records = tree.find("GuidDictionary"), tree.find("Records")
        if dictionary is None or records is None:
            raise ValueError("Missing saved Aspen ADS dictionary or records.")
        fields = {r.get("value"): r.get("int") for r in dictionary}
        self.fields = {key: fields.get(uid) for key, uid in _UUID.items()}
        self.records = {}
        self.by_type = defaultdict(list)
        for r in records:
            _unique(self.records, r.get("guid"), r)
            self.by_type[r.get("typeGuid")].append(r)

    def values(self, record):
        return {v.get("guid"): v.text for v in record.findall("Value")}

    def number(self, vals, field, optional=False):
        value = vals.get(self.fields[field])
        if optional and value is None:
            return None
        return _number(value, field)

    def entries(self, record, field):
        a = record.find(f"Array[@guid='{self.fields[field]}']")
        if a is None:
            raise ValueError(f"Missing Aspen {field} for record {record.get('guid')}.")
        return ((e.findtext("Key"), e.findtext("Val")) for e in a)

    def components(self, info):
        components, by_guid = {}, {}
        for r in self.by_type[self.fields["component_record"]]:
            v = self.values(r)
            cid = v[self.fields["component_id"]]
            if cid not in info:
                raise ValueError(f"Missing supported Aspen database identity for component {cid}.")
            _unique(components, cid, dict(id=cid, alias=v[self.fields["alias"]],
                                         molecular_weight_kg_kmol=self.number(v, "mw"), **info[cid]))
            by_guid[r.get("guid")] = cid
        if not components or set(components) != set(info):
            raise ValueError("Aspen component inventories differ between input and saved results.")
        # COMPONENTS MAIN determines the user-facing order in Aspen exports.
        # ADS object order can differ even within the same saved file.
        return [components[cid] for cid in info], by_guid

    def substreams(self, r, by_guid):
        rates = dict.fromkeys(by_guid.values(), 0.0)
        details = {}
        for substream, guid in self.entries(r, "substreams"):
            if substream not in {"MIXED", "CISOLID"}:
                raise ValueError(f"Unsupported Aspen substream {substream}.")
            sub = self.records[guid]
            values = self.values(sub)
            molar = self.number(values, "molar_flow") * 3600
            if molar == 0:
                continue
            fractions = {}
            for comp, value in self.entries(sub, "mole_fractions"):
                _unique(fractions, by_guid[comp], _number(value, "mole fraction"))
            if min(fractions.values(), default=-1) < 0 or not isclose(fsum(fractions.values()), 1., abs_tol=1e-8):
                raise ValueError(f"Incomplete mole fraction vector in saved substream {guid}.")
            subflows = {}
            for cid, fraction in fractions.items():
                rates[cid] += molar * fraction
                if fraction:
                    subflows[cid] = molar * fraction
            mass = self.number(values, "mass_flow")
            h = self.number(values, "h_j_kg", optional=True)
            s = self.number(values, "s_j_kg_k", optional=True)
            _unique(details, substream, dict(
                molar_kmol_h=molar, mass_kg_h=mass * 3600, component_flows_kmol_h=subflows,
                enthalpy_kw=mass * h / 1000 if h is not None else None,
                entropy_kw_k=mass * s / 1000 if s is not None else None,
            ))
        return rates, details

    def active_phases(self, vals, substreams):
        """Identify nonzero native phases; an incomplete phase inventory is unknown."""
        labels = []
        if "MIXED" in substreams:
            fluid_mass = 0.
            for field, label in (("vapor_record", "Vapor"), ("liquid1_record", "Liquid1"),
                                 ("liquid2_record", "Liquid2")):
                pointer = vals.get(self.fields[field])
                if pointer is None:
                    continue
                state = self.values(self.records[pointer])
                mass = self.number(state, "phase_mass_flow", optional=True)
                if mass and mass > 0:
                    fluid_mass += mass * 3600
                    labels.append(label)
            # The independently saved phase and MIXED masses can differ at
            # sub-ppm precision; the phase labels only describe the inventory.
            if not isclose(fluid_mass, substreams["MIXED"]["mass_kg_h"], rel_tol=1e-6, abs_tol=1e-10):
                return None
        if "CISOLID" in substreams:
            labels.append("CISOLID")
        return labels

    def streams(self, headers, by_guid):
        rows = {}
        for r in self.by_type[self.fields["material_record"]]:
            v = self.values(r)
            sid = v[self.fields["stream_id"]]
            composition, substreams = self.substreams(r, by_guid)
            mass = self.number(v, "mass_flow")
            molar = self.number(v, "molar_flow") * 3600
            if mass < 0 or molar < 0 or not isclose(fsum(composition.values()), molar, rel_tol=1e-8, abs_tol=1e-8):
                raise ValueError(f"Invalid bulk or component flows for {sid}.")
            if (mass == 0) != (molar == 0):
                raise ValueError(f"Inconsistent zero mass and molar flow for {sid}.")
            temp = self.number(v, "temperature", optional=mass == 0)
            pressure = self.number(v, "pressure", optional=mass == 0)
            _unique(rows, sid, dict(**headers[sid], temperature_k=temp, pressure_pa=pressure,
                                   enthalpy_kw=mass * self.number(v, "h_j_kg") / 1000 if mass else 0.,
                                   entropy_kw_k=mass * self.number(v, "s_j_kg_k") / 1000 if mass else 0.,
                                   mass_kg_h=mass * 3600, molar_kmol_h=molar,
                                   component_flows_kmol_h=composition,
                                   substreams=substreams, active_phases=self.active_phases(v, substreams),
                                   **{k: self.number(v, k, optional=True) for k in
                                      ("vapor_fraction", "liquid_fraction", "solid_fraction")}))
        if set(rows) != set(headers):
            raise ValueError("Aspen material stream inventories differ between ADS and text results.")
        return [rows[sid] for sid in headers]


def _transfer(name, kind, rate, **fields):
    return dict(id=name, kind=kind, rate_kw=rate, temperature_k=None, **fields)


def _block_results(ads, blocks, links):
    saved = {}
    for r in ads.records.values():
        bid = r.findtext(f"Value[@guid='{ads.fields['block_id']}']")
        if bid in blocks:
            _unique(saved, bid, ads.values(r))
    if set(saved) != set(blocks):
        raise ValueError("Missing saved Aspen equipment results.")
    for bid, b in blocks.items():
        ports = b.pop("energy_inlets") + b.pop("energy_outlets")
        kind = b["type"]
        v = saved[bid]
        transfers = b["transfers"]

        def kw(key):
            return ads.number(v, key) / 1000

        if kind in _HEAT_MODELS:
            net = kw("heat_w")
            for sid in ports:
                link = links[sid]
                if link["kind"] != "heat":
                    raise ValueError(f"Unsupported work connection at {bid}.")
                incoming = link["saved_rate_kw"] * (1 if link["to_block"] == bid else -1)
                # DSET links are rounded. Prefer exact ADS equipment duty when
                # this link alone supplies the entire duty.
                if len(ports) == 1 and isclose(net, incoming, rel_tol=6e-9, abs_tol=1e-10):
                    incoming = net
                internal = bool(link["from_block"] and link["to_block"])
                transfers.append(_transfer(sid, "heat", incoming, link_id=sid if internal else None))
            unlinked = net - fsum(t["rate_kw"] for t in transfers)
            if unlinked:
                transfers.append(_transfer("duty", "heat", unlinked))
        elif kind in {"PUMP", "COMPR"}:
            shaft = kw("pump_shaft_w" if kind == "PUMP" else "compressor_shaft_w")
            loss = shaft - kw("pump_fluid_w") if kind == "PUMP" else kw("compressor_loss_w")
            if len(ports) > 1 or any(links[sid]["kind"] != "work" for sid in ports):
                raise ValueError(f"Unsupported energy connection at {bid}.")
            name = ports[0] if ports else "shaft"
            link_id = None
            if ports:
                link = links[name]
                if not isclose(abs(shaft), abs(link["saved_rate_kw"]), rel_tol=6e-9, abs_tol=1e-8):
                    raise ValueError(f"Work stream {name} disagrees with {bid} shaft power.")
                if link["from_block"] and link["to_block"]:
                    link_id = name
            transfers.append(_transfer(name, "work", shaft, link_id=link_id))
            if loss:
                transfers.append(_transfer("drive_loss", "heat", -loss))
        elif ports:
            raise ValueError(f"Unsupported energy port on {kind} block {bid}.")
        elif kind == "HEATX":
            b["internal_heat_transfer_kw"] = kw("heat_w")
        elif kind == "RADFRAC":
            for service in ("condenser", "reboiler"):
                transfers.append(_transfer(service, "heat", kw(service + "_w"),
                                           process_temperature_k=ads.number(v, service + "_k", optional=True)))
        elif kind not in _PASSIVE_MODELS:
            raise ValueError(f"Unsupported Aspen block model {kind} at {bid}.")


def read_aspen(path):
    """Return the normalized stream, block, component and energy-link inventory.

    The signed transfer rates are positive into the equipment. The separate
    energy_links table preserves Aspen's own saved link value and port direction.
    An unmeasured heat boundary temperature is returned as None.
    """
    path = Path(path)
    source = path.read_bytes()
    member, content = _backup(path, source)
    if not re.match(rb'^MM "41\.0"(?:\s|$)', content):
        raise ValueError("Unsupported Aspen Plus backup version. This reader supports 41.0.")
    text, marker, rest = content.partition(b"$_ADS_FILE")
    payload, end, _ = rest.partition(b"$_END_ADS_FILE")
    if not marker or not end:
        raise ValueError("Missing embedded Aspen ADS results. Save a solved case including results.")
    try:
        tree = ET.fromstring(zlib.decompress(base64.b64decode(payload), -15))
    except (ValueError, zlib.error, ET.ParseError) as error:
        raise ValueError("Invalid embedded Aspen ADS XML results.") from error
    text = text.decode("latin-1").replace("\r\n", "\n")
    flat = " ".join(text.split(" DSET ", 1)[0].split())
    blocks = _topology(flat)
    headers, links, warnings = _result_headers(text)
    _check_topology(blocks, headers, links)
    ads = _ADS(tree)
    try:
        comps, by_guid = ads.components(_component_info(flat))
        streams = ads.streams(headers, by_guid)
        _block_results(ads, blocks, links)
    except KeyError as error:
        raise ValueError(f"Missing or inconsistent Aspen record {error}.") from error
    sha = hashlib.sha256(content).hexdigest()
    section_options = _reference_option_hints(flat, streams)
    metadata = dict(source_filename=path.name, source_member=member, backup_version="41.0",
                    backup_sha256=sha, source_sha256=hashlib.sha256(source).hexdigest(),
                    source_properties="Aspen saved ADS SI bulk properties", section_property_options=section_options,
                    enthalpy_entropy_basis="aspen_formation" if sha in
                    {_AUDITED_FORMATION_BACKUP, EXERGY_REPORT_BACKUP} else "unverified")
    return dict(metadata=metadata, components=comps, streams=streams, blocks=list(blocks.values()),
                energy_links=list(links.values()), warnings=warnings)


def assert_cache_units(model):
    """Reject recognizable display-unit block caches; never guess a conversion.

    Some Aspen COM backup exports mix SI material records with displayed MW/C
    equipment records under the same property UUIDs. Those fields carry no
    reliable unit label in ADS. Require an explicit COM/report normalization
    when independent stream states demonstrate that inconsistency.
    """
    streams = {s["id"]: s for s in model["streams"]}
    columns, exchangers = [], []
    for block in model["blocks"]:
        if block["type"] == "RADFRAC":
            temperatures = [streams[s].get("temperature_k") for s in block["outlets"]]
            temperatures = [t for t in temperatures if t is not None]
            shifted = 0
            for transfer in block["transfers"]:
                t = transfer.get("process_temperature_k")
                if t is None or transfer["id"] not in {"condenser", "reboiler"}:
                    continue
                matches = lambda candidate: any(isclose(candidate, value, rel_tol=1e-6, abs_tol=1e-3)
                                                for value in temperatures)
                if not matches(t) and matches(t + 273.15):
                    shifted += 1
            if shifted == 2:
                columns.append(block["id"])
        elif block["type"] == "HEATX":
            duty = block.get("internal_heat_transfer_kw")
            if duty is None or abs(duty) < 1e-12:
                continue
            differences = [abs(streams[o]["enthalpy_kw"] - streams[i]["enthalpy_kw"])
                           for i in block["inlets"] for o in block["outlets"]
                           if streams[i].get("enthalpy_kw") is not None
                           and streams[o].get("enthalpy_kw") is not None]
            matches = lambda candidate: any(isclose(abs(candidate), value, rel_tol=1e-5, abs_tol=1e-3)
                                            for value in differences)
            if not matches(duty) and matches(duty * 1e6):
                exchangers.append(block["id"])
    if columns or len(exchangers) >= 2:
        affected = ", ".join(columns + exchangers)
        raise ValueError(
            "Aspen block cache contains inconsistent display units relative to SI stream states "
            f"({affected}). Normalize using explicit native COM/report units and load the "
            "validated normalized JSON; no automatic unit conversion was applied."
        )
