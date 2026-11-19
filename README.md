# Aspen Energy Exergy

A Python library by **Mohammadreza Jahanbakhsh** for steady-state energy and exergy analysis of saved Aspen Plus results. It calculates stream enthalpy and physical, chemical and total exergy, equipment and subsystem balances, signed exergy destruction, entropy generation and selected product performance. Every stream records the source and quality of its exergy split.

Version 0.5.1 supports the built-in balance engine and an imported-state engine using **ExerPy 0.2.2** and **TESPy 0.12.0.post1**. The repository includes a small hydrogen-heater example and automated tests. Reading supported saved results and running the example requires Python 3.10 or later. Generating a new Aspen simulation requires your Aspen installation and licenses.

## Install and run

Clone the repository, install the package into your Python environment and analyze the included example.

```sh
git clone https://github.com/mrezajbh/aspen-energy-exergy.git
cd aspen-energy-exergy
python -m pip install .
python -m aspen_energy_exergy examples/heater.json -o heater-results.zip
```

The example uses an ideal-gas hydrogen heater with constant heat capacity. Its first-law residual is zero within numerical precision, and its signed exergy destruction is approximately **0.122 kW**. Default physical and chemical splits use external property estimates and carry the `estimated` quality label.

```python
from aspen_energy_exergy import analyze, write_results

report = analyze("examples/heater.json")
print(report["plant"]["energy_residual_kw"])
print(report["plant"]["exergy_destruction_kw"])
print(report["streams"]["PRODUCT"]["split_quality"])
write_results(report, "heater-results.zip")
```

The ZIP contains a readable summary, `analysis.json` and CSV tables for streams, equipment, heat/work transfers, the plant and quality findings. Configured groups and fuel selections add subsystem and fuel-performance tables. Rates use kW, entropy rates use kW/K, flows use kg/h or kmol/h, temperatures use K and pressures use Pa. A blank CSV field or JSON `null` indicates an unavailable quantity. Review `quality.csv` for missing inputs, discrepancies, phase limitations and signed negative destruction.

## Use ExerPy and TESPy balances

The `exerpy-tespy` engine executes the bundled upstream balance routines with the supplied stream properties and declared heat/work boundaries. Select it together with `estimate_splits=False` in Python or `--native-splits-only` at the command line.

```sh
python -m aspen_energy_exergy examples/heater.json --native-splits-only --balance-engine exerpy-tespy -o heater-external-results.zip
```

```python
report = analyze(
    "examples/heater.json",
    estimate_splits=False,
    balance_engine="exerpy-tespy",
)
print(report["plant"]["energy_residual_kw_engine"])
print(report["plant"]["exergy_destruction_kw_engine"])
```

This example provides the properties needed for total-exergy and energy balances. Its physical and chemical splits remain unavailable in this mode because the example supplies no native EXERGYFL or ambient flash. Supply those results to obtain a complete native split.

The bundled ExerPy and TESPy sources retain their licenses, source hashes and patches. The patches make unused property-backend imports optional. This engine evaluates imported-state balances, and a runtime import guard checks that CoolProp is absent. Use a fresh Python process when selecting it. The `external_engines` result records upstream versions, methods and comparisons with the built-in calculations. Agreement tests the balance implementation. The input simulation's thermodynamic consistency still requires review.

## Analyze an Aspen flowsheet

Replace the example path with your solved `.apwz` or `.bkp` file.

```sh
python -m aspen_energy_exergy process.apwz --settings settings.json -o process-results.zip
```

The native reader supports Aspen Plus 41.0 archives containing compatible embedded ADS results and MIXED/CISOLID substreams. Recognized equipment includes heaters, two-sided HEATX units, pumps, compressors, valves, mixers, splitters, separators, FLASH2, RADFRAC and selected reactors. Unsupported versions, units, substreams and layouts produce an explicit error. For other Aspen versions, map a solved export to the normalized structure in `examples/heater.json`. `load_case(path)` returns that same structure as a Python dictionary.

A normalized case has `metadata`, `components`, `streams` and `blocks`.

- Components need an ID and elemental composition, such as `{"C": 1, "H": 4}`. Supply a CAS number or unambiguous chemical name for external property calculations.
- Streams identify source and destination blocks, flow rates, operating properties, component rates and the Aspen property method. External ends use JSON `null`. Enthalpy and entropy are whole-material rates including active substreams. Vapor fraction can disambiguate saturated water.
- Blocks list material inlets and outlets and individual heat/work transfers. Transfers are positive into the block. A shared energy stream uses the same `link_id` on its two transfers, with opposite rates. A two-sided exchanger can record `internal_heat_transfer_kw`, which cancels inside that block's boundary.
- The elemental-reference reconstruction requires Aspen's elemental formation convention for both enthalpy and entropy. Set `metadata.enthalpy_entropy_basis` to `aspen_formation` only after confirming that convention. An unverified datum leaves the affected computed quantities unavailable. A settings file can declare a checked native formation basis with `confirm_formation_datum: true`.

### Import native physical exergy

A matching Aspen material report workbook can supply operating states and EXERGYFL. The library checks it against the associated archive. A matching text report adds a second exported-data comparison.

```python
report = analyze(
    "process.apwz",
    native_exergy="material-report.xlsx",
    report_text="process-report.txt",
    settings="settings.json",
    estimate_splits=False,
    balance_engine="exerpy-tespy",
)
```

The workbook importer expects the supported Aspen material-report layout. It records source locations, hashes, substream scope and warnings in the result. Aspen's exported reference settings and property conventions must be compatible with the analysis environment.

For a compact CSV exchange, generate a case-specific return sheet and fill the available EXERGYFL values in kW.

```sh
python -m aspen_energy_exergy process.apwz --write-exergy-template Aspen_EXERGYFL_return.csv
python -m aspen_energy_exergy process.apwz --native-exergy Aspen_EXERGYFL_return.csv --native-splits-only -o checked-results.zip
```

Enter the report location, confirm the operating state, and supply the actual reference temperature and pressure. Choose `whole_stream`, `substreams` or `mixed_with_saved_graphite` in the `scope` field and complete the corresponding value columns. Leave outstanding rows empty. A numerical zero counts as a returned value. The importer preserves provenance and checks conflicting direct rates against any matching ambient flash.

### Request ambient reference flashes

The planner groups outstanding streams by compatible composition, substream allocation and saved property settings. Each group identifies a representative recipe and possible peers. Confirm the complete settings in Aspen before reusing a flash for those peers.

```sh
python -m aspen_energy_exergy process.apwz --write-flash-plan ambient-requests.zip
python -m aspen_energy_exergy process.apwz --references flash_returns.csv --native-splits-only -o checked-results.zip
python -m aspen_energy_exergy process.apwz --references flash_returns.csv --write-flash-plan remaining-requests.zip
```

The request ZIP contains instructions, compositions, groups, a CSV return sheet and a fillable JSON alternative. Flash the representative at **298.15 K and 101325 Pa** with matching composition, property method, interaction data, substream allocation and allowed phases, including free-water settings. In `flash_returns.csv`, enter whole-material molar enthalpy in J/mol and entropy in J/(mol K), actual reference conditions, the flash ID and its report location. Aspen kJ/kmol and kJ/(kmol K) have the same numerical values. Confirm the representative and retain only verified peers in `applies_to`. Return fields for outstanding flashes can stay blank. Repeat `--references` to combine nonoverlapping batches.

The equivalent Python API provides `plan_reference_flashes` and `write_flash_plan`. A `.json` destination writes the fillable JSON alone. For independently prepared reference JSON, supply a `case_sha256` matching the case and a `streams` array. Each record needs its operating stream `id`, `aspen_reference_id`, `source_location`, `T_K`, `P_Pa`, `property_method`, `same_aspen_property_options: true`, `phase_convention: "equilibrium"`, ordered `component_ids` and `component_flows_kmol_h`. Choose one quantity basis.

| `quantity_basis` | Enthalpy field | Entropy field |
| --- | --- | --- |
| `specific_SI` | `h0_J_kg` | `s0_J_kgK` |
| `molar_SI` | `h0_J_mol` | `s0_J_molK` |
| `flow_SI` | `H0_kW` | `S0_kW_K` |

### Set heat boundaries and products

Pass a settings dictionary to `analyze` or a JSON file to `--settings`. Each nonzero heat transfer needs a thermodynamic boundary temperature to calculate its exergy. Every override records the transfer and the evidence or assumption behind it.

```python
settings = {
    "name": "Hydrogen heater with measured steam temperature",
    "transfers": [{
        "block": "HEATER",
        "id": "steam",
        "temperature_k": 420.0,
        "basis": "Steam header measurement",
    }],
    "product_streams": ["PRODUCT"],
    "groups": {"heater_package": ["HEATER"]},
}
report = analyze("examples/heater.json", settings=settings)
print(report["groups"]["heater_package"]["exergy_destruction_kw"])
```

Use a declared `kind` override when a process uses an Aspen heat port to represent electrical work. For a sensible heat contact with constant heat-capacity rate, `temperature_range_k` accepts two positive endpoint temperatures in place of `temperature_k`. The calculation integrates `dQ/T` using their logarithmic mean. Intervals crossing the environmental temperature record incoming and outgoing exergy separately. Model phase changes with separate sensible and latent transfers and their applicable temperatures.

The bundled `biogas-ft` settings describe a particular research case and select automatically only when its backup hash matches. Other native files use their saved transfer inventory by default. Inspect every chosen profile's boundary assumptions when adapting a case. Set `settings="as-saved"` to use the raw inventory. The profile source records each assumption in `src/aspen_energy_exergy/profiles/biogas-ft.json`.

`groups` maps subsystem names to block IDs. Each result uses that group's external boundary. Overlapping groups are independent views. Plant totals cancel internal streams and transfers directly. `product_streams` selects external products for the reported share of incoming exergy.

For fuel-product performance, configure `fuel_feeds`, `product_streams` and any supply/return `utility_pairs`. The efficiency divides product lower-heating-value flow by declared feed LHV, gross work input, external heat input and consumed steam energy. The result records this boundary. Pure-component standard combustion data come from `thermo` and `chemicals`. [ChEDL chemical constants](https://thermo.readthedocs.io/thermo.chemical_package.html)

For multiple scenarios, load the case once with `model = load_case(path)` and pass it to each analysis. Each call applies its own settings and native returns. Property estimates at identical thermodynamic inputs are cached within the Python process.

## Thermodynamic basis and limitations

The environmental reference is **298.15 K and 101325 Pa**. With formation-referenced Aspen enthalpy and entropy, total stream exergy in kW is `H − T₀S + Σ(nₑ b°ₑ)/3600`, where `nₑ` is in kmol of atoms/h and `b°ₑ` is in J/mol of atoms. The boundary includes thermomechanical and chemical exergy. Bulk kinetic and gravitational energy are outside its scope. [Tang et al., equation 5 and table 1](https://www.aidic.it/cet/21/83/089.pdf)

A matching ambient flash gives physical exergy as `(H − H₀) − T₀(S − S₀)`. Native EXERGYFL provides a direct physical rate. Chemical exergy is the residual of the computed total after subtracting the accepted physical rate. Stream provenance identifies the supplied source and the accepted substream scope.

In the default engine, missing native splits can use external property estimates. Water and steam use IAPWS-IF97. Fluid mixtures use ChEDL property routines with matched operating/reference calculations. Compatible streams can share an estimated ambient molar chemical exergy, with independent estimates retained for comparison. Failed mixture flashes may use an explicitly labelled Wilson K-value approximation. Review `split_method`, `split_quality` and `notes` for phase assumptions and extrapolation. Disable estimation with `estimate_splits=False` or `--native-splits-only`. [ChEDL phases and flashes](https://thermo.readthedocs.io/tutorial_phases_and_flash.html) · [IAPWS industrial formulation](https://iapws.org/technical-guidance/release/IF97-Rev)

For each boundary, the first-law residual is material enthalpy in minus material enthalpy out plus signed heat and work. Each heat contact contributes `Q(1 − T₀/Tb)` to exergy and `Q/Tb` to entropy transfer. Shaft and electrical work carry their full rate as exergy. Energy closure, mass and elemental residuals, entropy generation and signed destruction appear together. Negative destruction flags a boundary or model discrepancy. Convergence and complete exported data alone do not establish thermodynamic consistency, especially across property-method interfaces or conceptual separations.

## Tests

```sh
python -m unittest discover -s tests -v
```

The portable tests cover the public API, balances, properties, native-return validation, flash planning and ExerPy/TESPy integration. Additional regression tests detect optional local Aspen study fixtures and skip when those files are unavailable. The included heater example supports a complete installation and command-line check.

## Author and license

**Mohammadreza Jahanbakhsh** maintains this library. The project is free software distributed under **GNU General Public License version 3 only (`GPL-3.0-only`)**. You may use, study, modify and redistribute it under those terms. Distribution of covered modifications requires the corresponding source and license notices. See [LICENSE](LICENSE).

Third-party libraries retain their authorship and license terms. The bundled ExerPy and TESPy sources include their original notices and modification records. Dependency and attribution details are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
