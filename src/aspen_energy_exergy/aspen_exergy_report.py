"""Read and reconcile an Aspen Plus v41 workbook, printed report and saved case.

The source workbook provides full-precision extensive physical exergy and the
operating state of each stream, including MIXED and conventional solid flows.
Only the Python standard library is needed. No property model is evaluated.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import posixpath
import re
import xml.etree.ElementTree as ET
from zipfile import ZipFile

from .aspen import read_aspen

CAL_PER_S_TO_KW = 0.0041868  # Aspen uses the international-table calorie
SHEET = "Aspen Data Tables"

def _workbook_rows(path):
    with ZipFile(path) as z:
        catalog = ET.fromstring(z.read("xl/workbook.xml"))
        sheet = next(s for s in catalog.findall(".//{*}sheet") if s.get("name") == SHEET)
        rid = sheet.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = next(r.get("Target") for r in rels if r.get("Id") == rid)
        name = target.lstrip("/") if target.startswith("/") else posixpath.join("xl", target)
        strings = ["".join(s.itertext()) for s in ET.fromstring(z.read("xl/sharedStrings.xml"))]
        rows = {}
        with z.open(name) as xml:
            for _, row in ET.iterparse(xml, events=("end",)):
                if not row.tag.endswith("}row"):
                    continue
                cells = {}
                for cell in row:
                    col = cell.get("r", "").rstrip("0123456789")
                    raw = cell.findtext("{*}v")
                    if raw is not None:
                        cells[col] = strings[int(raw)] if cell.get("t") == "s" else raw
                rows[int(row.get("r"))] = cells
                row.clear()
    return rows


def _find_row(rows, label, start=0, end=10**9):
    return next(n for n, r in rows.items() if start <= n < end and r.get("C") == label)


def _section(rows, name, end):
    start = _find_row(rows, name)
    layout = {}
    for n in range(start + 1, end):
        row = rows[n]
        if "C" in row:
            layout.setdefault(row["C"], n)
    expected = {"Temperature":"C", "Pressure":"bar", "Mass Flows":"kg/hr",
                "Mole Flows":"kmol/hr", "Enthalpy Flow":"cal/sec",
                "Mass Entropy":"cal/gm-K", "Molar Entropy":"cal/mol-K",
                "Exergy flow rate":"cal/sec"}
    for label, unit in expected.items():
        if rows[layout[label]].get("D") != unit:
            raise ValueError(f"Unsupported units in {name} / {label}")
    layout["component_rows"] = {
        rows[n]["C"]:n for n in range(layout["Mole Flows"]+1,layout["Mole Fractions"])
        if rows[n].get("D") == "kmol/hr"
    }
    return layout


def _number(rows, n, col):
    raw = rows.get(n, {}).get(col)
    if raw is None:
        return None
    value = float(raw)
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite input at {col}{n}")
    return value


def _state(rows, layout, col, zero):
    def value(key):
        return _number(rows, layout.get(key), col)

    mass, moles = value("Mass Flows"), value("Mole Flows")
    rates = {cid: _number(rows, n, col) for cid, n in layout["component_rows"].items()}
    h, smass, smolar = value("Enthalpy Flow"), value("Mass Entropy"), value("Molar Entropy")
    if zero:
        if any(x not in (None, 0.0) for x in (mass,moles,h,*rates.values())):
            raise ValueError(f"Nonzero extensive value for zero-flow column {col}")
        mass, moles, h, entropy = 0., 0., 0., 0.
    else:
        if None in (mass,moles,h,smass,smolar) or any(x is None or x < 0 for x in rates.values()):
            raise ValueError(f"Incomplete flowing operating state at {col}{layout['Mole Flows']}")
        if not math.isclose(math.fsum(rates.values()),moles,rel_tol=1e-11,abs_tol=1e-10):
            raise ValueError(f"Species do not sum to reported molar flow in {col}")
        entropy = mass * smass * 4.1868 / 3600
        if not math.isclose(entropy, moles*smolar*4.1868/3600,rel_tol=1e-10,abs_tol=1e-10):
            raise ValueError(f"Molar and mass entropy rates differ in {col}")
    t, p = value("Temperature"), value("Pressure")
    return {
        "enthalpy_kw": h * CAL_PER_S_TO_KW, "entropy_kw_k": entropy,
        "temperature_k": t + 273.15 if t is not None else None,
        "pressure_pa": p * 1e5 if p is not None else None,
        "mass_kg_h": mass, "molar_kmol_h": moles,
        "component_flows_kmol_h": {cid:x for cid,x in rates.items() if x},
        **{kind.lower()+"_fraction": value("Molar "+kind+" Fraction")
           for kind in ("Vapor", "Liquid", "Solid")},
    }


def _text_mixed(path):
    result = {}
    ids = []
    substream = None
    with Path(path).open(encoding="latin1") as report:
        for line_no, line in enumerate(report, 1):
            if line.startswith(" STREAM ID "):
                ids = line.split()[2:]
                substream = None
            section = re.match(r"^ (\S+)\s+SUBSTREAM PROPERTIES:", line)
            if section:
                substream = section[1]
            if not line.startswith(" EXERGYFL "):
                continue
            if substream != "MIXED":
                continue
            _, unit, *values = line.split()
            if unit.upper() != "KW" or len(values) != len(ids):
                raise ValueError(f"Unrecognized report units or stream columns at line {line_no}")
            for sid, token in zip(ids,values):
                if sid in result:
                    raise ValueError(f"Duplicate printed MIXED EXERGYFL for {sid}")
                if token.upper() == "MISSING":
                    result[sid] = (None,0.,line_no)
                else:
                    # Aspen prints compact exponents, e.g. 6.2753-02 or 1.5966+04.
                    token = re.sub(r"(?<=\d)([+-]\d+)$", r"E\1", token)
                    mant, _, exp = token.upper().partition("E")
                    places = len(mant.partition(".")[2])
                    half_step = 0.5 * 10**(int(exp or 0)-places)
                    result[sid] = (float(token),half_step,line_no)
    return result


def parse_export(workbook, text_report, case_path, *, case_model=None):
    """Return source-traceable native physical exergy, composition and H/S by stream."""
    workbook, case_path = Path(workbook), Path(case_path)
    text_report = Path(text_report) if text_report is not None else None
    sources = [workbook, case_path] + ([text_report] if text_report else [])
    source_hashes = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    model = case_model if case_model is not None else read_aspen(case_path)
    if model["metadata"]["source_sha256"] != source_hashes[case_path.name]:
        raise ValueError("Previously read model differs from supplied Aspen case")
    rows = _workbook_rows(workbook)
    header = _find_row(rows,"Stream Name")
    columns = {c:i for c,i in rows[header].items() if c not in ("C","D")}
    if len(set(columns.values())) != len(columns):
        raise ValueError("Duplicate workbook stream identifier")
    mixed_kw = _find_row(rows,"(MIXED) Exergy flow rate")
    solid_kw = _find_row(rows,"(SOLIDS) Exergy flow rate")
    if any(rows[n].get("D") != "kW" for n in (mixed_kw,solid_kw)):
        raise ValueError("Expected EXERGYFL substream rates in kW")
    total = _section(rows,"Total Stream",_find_row(rows,"MIXED Substream"))
    mixed = _section(rows,"MIXED Substream",mixed_kw+1)
    solid = _section(rows,"CISOLID Substream",solid_kw+1)
    case = {r["id"]:r for r in model["streams"]}
    species = {r["id"] for r in model["components"]}
    from_row, to_row = _find_row(rows,"From"), _find_row(rows,"To")
    printed = _text_mixed(text_report) if text_report else None
    if set(columns.values()) != set(case) or (printed is not None and set(case) != set(printed)):
        raise ValueError("Workbook, report and archive stream inventories differ")
    if any(set(layout["component_rows"]) != species for layout in (total,mixed,solid)):
        raise ValueError("Workbook and case component inventories differ")
    streams = {}
    maxima = {}
    text_max = 0.
    excess_max = 0.
    closure_max = 0.
    def check(field,sid,a,b,atol,rtol):
        if a is None and b is None:
            return
        if a is None or b is None or not math.isclose(a,b,abs_tol=atol,rel_tol=rtol):
            raise ValueError(f"Case/workbook {field} differs for {sid}: {a} vs {b}")
        delta = abs(a-b)
        if field not in maxima or delta > maxima[field]["absolute_difference"]:
            maxima[field] = {"stream":sid, "absolute_difference":delta}
    for col,sid in columns.items():
        for n, field in ((from_row,"from_block"),(to_row,"to_block")):
            if rows[n].get(col) != case[sid][field]:
                raise ValueError(f"Workbook and case {field} differ for {sid}")
        subs = {}
        for label,layout in (("MIXED",mixed),("CISOLID",solid)):
            n = _number(rows,layout["Mole Flows"],col)
            if n is None:
                raise ValueError(f"Missing substream molar flow for {sid}/{label}")
            subs[label] = _state(rows,layout,col,n == 0)
        zero = not any(s["molar_kmol_h"] for s in subs.values())
        op = _state(rows,total,col,zero)
        for field in ("enthalpy_kw","entropy_kw_k","mass_kg_h","molar_kmol_h"):
            sub_total = math.fsum(s[field] for s in subs.values())
            if not math.isclose(op[field],sub_total,abs_tol=1e-8,rel_tol=1e-11):
                raise ValueError(f"Substream sum disagrees with workbook whole {field} for {sid}")
        for field,atol,rtol in (("temperature_k",0.005,0.),("pressure_pa",0.01,0.),
                                ("enthalpy_kw",1e-5,1e-5), ("entropy_kw_k",1e-8,1e-5),
                                ("mass_kg_h",1e-6,1e-5), ("molar_kmol_h",1e-8,1e-5)):
            check(field,sid,op[field],case[sid][field],atol,rtol)
        for cid in species:
            flow = op["component_flows_kmol_h"].get(cid,0.)
            sub_total = sum(s["component_flows_kmol_h"].get(cid,0.) for s in subs.values())
            if not math.isclose(flow,sub_total,abs_tol=1e-10,rel_tol=1e-11):
                raise ValueError(f"Substream amounts differ from whole for {sid}/{cid}")
            check("component_kmol_h",sid+"/"+cid,flow,
                  case[sid]["component_flows_kmol_h"].get(cid,0.),0.0002,1e-5)
        op["substreams"] = subs
        mx, so = _number(rows,mixed_kw,col), _number(rows,solid_kw,col)
        for label,rate in (("MIXED",mx),("CISOLID",so)):
            if subs[label]["molar_kmol_h"] > 0 and rate is None:
                raise ValueError(f"Missing flowing substream EXERGYFL for {sid}/{label}")
            if subs[label]["molar_kmol_h"] == 0 and rate not in (0.,None):
                raise ValueError(f"Nonzero EXERGYFL for inactive {sid}/{label}")
        native = _number(rows,total["Exergy flow rate"],col)
        if native is None and not zero:
            raise ValueError(f"Missing flowing whole-stream exergy at {sid}")
        whole = (native or 0.) * CAL_PER_S_TO_KW
        closure = abs(whole-math.fsum(x for x in (mx,so) if x is not None))
        if closure > max(1e-9,abs(whole)*1e-12):
            raise ValueError(f"Whole-stream EXERGYFL disagrees with substreams for {sid}")
        closure_max = max(closure_max,closure)
        text_val,precision,text_line = printed[sid] if printed is not None else (None,None,None)
        if printed is not None and (mx is None) != (text_val is None):
            raise ValueError(f"Printed-report and workbook MIXED availability differs at {sid}")
        if text_val is not None:
            delta = abs(mx-text_val)
            if delta > precision + max(0.0001,abs(mx)*6e-6):
                raise ValueError(f"Printed-report and workbook MIXED values disagree at {sid}")
            text_max = max(text_max,delta)
            excess_max = max(excess_max,delta-precision)
        location = (f"{SHEET}!{col}{total['Exergy flow rate']},"
                    f"{col}{mixed_kw},{col}{solid_kw}")
        streams[sid] = {"mixed_kw":mx, "solid_kw":so, "whole_kw":whole,
                        "source_location":location, "source_column":col,
                        "native_whole_cal_s":native, "zero_flow":zero,
                        "text_mixed_kw":text_val, "text_report_line":text_line, "operating":op}
    summary = {
        "sha256":source_hashes, "text_report_checked":printed is not None,
        "topology_checked":True,
        "stream_count":len(streams),
        "flowing_streams":sum(not v["zero_flow"] for v in streams.values()),
        "mixed_reported":sum(v["mixed_kw"] is not None for v in streams.values()),
        "solid_reported":sum(v["solid_kw"] is not None for v in streams.values()),
        "no_flow_ids":[s for s,v in streams.items() if v["zero_flow"]],
        "negative_whole_kw":{s:v["whole_kw"] for s,v in streams.items() if v["whole_kw"]<0},
        "source_sheet":SHEET, "id_row":header, "whole_cal_s_row":total["Exergy flow rate"],
        "mixed_kw_row":mixed_kw, "cisolid_kw_row":solid_kw,
        "cal_s_to_kw":CAL_PER_S_TO_KW, "whole_minus_substreams_max_abs_kw":closure_max,
        "case_workbook_max_abs_differences":maxima, "component_count":len(species),
        "max_text_mixed_abs_difference_kw": text_max if printed is not None else None,
        "max_text_mixed_difference_beyond_display_rounding_kw":excess_max if printed is not None else None,
        "reference_temperature_k":298.15, "reference_pressure_pa":101325.,
        "reference_basis":"Aspen EXERGYFL documented default, consistent with the requested reference. The exports do not print an explicit environmental setting.",
        "reference_source":"https://doi.org/10.3390/en15010205",
        "state_basis":"Full-precision Aspen workbook for all stream states and physical exergy; archive and supplied printed report serve as independent checks.",
    }
    return {"streams":streams, "summary":summary}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("workbook","text_report","case_path","output_json"):
        parser.add_argument(arg,type=Path,**({"nargs":"?"} if arg=="text_report" else {}))
    args = parser.parse_args()
    if args.output_json.resolve() in {p.resolve() for p in (args.workbook,args.text_report,args.case_path) if p is not None}:
        parser.error("Output must differ from each source file")
    result = parse_export(args.workbook,args.text_report,args.case_path)
    args.output_json.write_text(json.dumps(result,indent=2,allow_nan=False),encoding="utf-8")
