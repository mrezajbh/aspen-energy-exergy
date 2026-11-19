"""Heat utilization and molecular-sieve design for the biogas-to-jet case."""

from math import fsum

from .balances import _transfer


def add_case_details(report):
    """Add case-specific reporting views to the already calculated balances."""
    streams = report["streams"]
    groups = report["groups"]
    heat = {(t["block"], t["id"]): t for t in report["transfers"] if t["kind"] == "heat"}
    hru = heat["HRU", "duty"]
    circuits = []
    for inlet, outlet, role in (
        ("DH-IN", "DH-OUT", "district"), ("DH-IN2", "DH-OUT2", "district"),
        ("DH-IN3", "DH-OUT3", "district"), ("DH-IN4", "DH-OUT4", "district"),
        ("DH-IN5", "DH-OUT5", "district"), ("DH-IN8", "DH-OUT8", "district"),
        ("DH-IN18", "DH-OUT18", "district"), ("CWS-01", "CWR-01", "cooling"),
        ("CWS-02", "CWR-02", "cooling"),
    ):
        a, b = streams[inlet], streams[outlet]
        block = a["to_block"]
        if not block or block != b["from_block"]:
            raise ValueError(f"Recipient water pair {inlet}/{outlet} spans different blocks")
        circuits.append(dict(block=block, role=role, inlet=inlet, outlet=outlet,
                             inlet_c=a["temperature_k"] - 273.15, outlet_c=b["temperature_k"] - 273.15,
                             heat_gain_kw=b["enthalpy_kw"] - a["enthalpy_kw"],
                             exergy_gain_kw=b["total_exergy_kw"] - a["total_exergy_kw"]))

    def credits(members):
        roles = [c for c in circuits if c["block"] in members]
        district = [c for c in roles if c["role"] == "district"]
        district_heat = fsum(c["heat_gain_kw"] for c in district)
        district_exergy = fsum(c["exergy_gain_kw"] for c in district)
        hru_heat = -hru["rate_kw"] if "HRU" in members else 0.
        hru_exergy = -hru["exergy_kw"] if "HRU" in members else 0.
        return dict(district_heat_kw=district_heat, district_heat_exergy_kw=district_exergy,
                    cooling_water_gain_kw=fsum(c["heat_gain_kw"] for c in roles if c["role"] == "cooling"),
                    hru_heat_kw=hru_heat, hru_heat_exergy_kw=hru_exergy,
                    recovered_heat_kw=district_heat + hru_heat,
                    recovered_heat_exergy_kw=district_exergy + hru_exergy)

    for group in groups.values():
        group.update(credits(group["members"]))
    plant_credits = credits(set(report["plant"]["members"]))
    report["heat_recovery"] = dict(
        water_circuits=circuits,
        hru_heat_kw=plant_credits["hru_heat_kw"], hru_heat_exergy_kw=plant_credits["hru_heat_exergy_kw"],
        hru_delivery_c=hru["temperature_k"] - 273.15,
        hru_inlet_c=streams["EFF-OUT"]["temperature_k"] - 273.15,
        hru_outlet_c=streams["COLED"]["temperature_k"] - 273.15,
        cold_sep_cooling_kw=-heat["COLD-SEP", "duty"]["rate_kw"],
        plant_district_heat_kw=plant_credits["district_heat_kw"],
        plant_district_heat_exergy_kw=plant_credits["district_heat_exergy_kw"],
        plant_recovered_heat_kw=plant_credits["recovered_heat_kw"],
        plant_recovered_heat_exergy_kw=plant_credits["recovered_heat_exergy_kw"],
        basis=("District heat and exergy are recipient water outlet minus inlet rates. "
               "HRU credit assumes full use at the specified uniform delivery temperature. "
               "The balances already include the water streams and HRU duty. Cooling-water gains "
               "and COLD-SEP cooling are reported separately."),
    )

    def water_ppmv(sid):
        stream = streams[sid]
        return stream["component_flows_kmol_h"].get("H2O", 0.) / stream["molar_kmol_h"] * 1e6

    water_mw = next(c["molecular_weight_kg_kmol"] for c in report["components"] if c["id"] == "H2O")
    removed = streams["WAT-D"]["component_flows_kmol_h"]["H2O"] * water_mw
    installed, capacity, adsorption_h = 5540., 0.10, 8.
    required = removed * adsorption_h / capacity
    outlet = water_ppmv("DRYED")
    regeneration = _transfer("DRYER", dict(id="design_regeneration", kind="heat", rate_kw=120.,
                             temperature_k=493.15), report["environment"]["temperature_k"])
    report["dryer_design"] = dict(
        outlet_ppmv=outlet, outlet_limit_ppmv=5., outlet_meets_limit=outlet <= 5.,
        hx04_inlet_bulk_ppmv=water_ppmv("TO-CHILL"), water_removed_kg_h=removed,
        bed_count=2, installed_mass_per_bed_kg=installed, required_mass_per_bed_kg=required,
        bed_mass_shortfall_kg=max(0., required - installed), working_capacity_kg_per_kg=capacity,
        adsorption_h=adsorption_h, regeneration_h=8., purge_co2_kmol_h=10.,
        regeneration_heat_kw=regeneration["rate_kw"], regeneration_temperature_c=220.,
        regeneration_heat_exergy_kw=regeneration["exergy_kw"],
        ammonia_inlet_c=streams["AMMONIA"]["temperature_k"] - 273.15,
        user_frost_c_at_18_bar=-43., user_frost_c_at_1_atm=-65.,
        basis=("DRYED and TO-CHILL moisture use bulk wet mole fractions. WAT-D gives the separated "
               "water load. The user supplied the alternating 3A/4A beds, cycle duration, capacity, "
               "CO2 purge, frost points and regeneration heat. The Aspen package models continuous "
               "separation. Cyclic regeneration heat and purge are separate design requirements. "
               "The export supplies no regeneration purge outlet or spent-media state."),
    )
