# Licensing and third-party notices

Aspen Energy Exergy is maintained by Mohammadreza Jahanbakhsh. Copyright (c) 2026 Mohammadreza Jahanbakhsh applies to the original project code, documentation, and modifications. The project is distributed under the GNU General Public License, version 3 only (`GPL-3.0-only`). The complete terms are in [LICENSE](LICENSE).

The required `iapws` dependency uses GPLv3. Distributing the integrated project under GPLv3 accommodates that dependency and the compatible MIT and BSD licenses used by the other libraries. This license permits use, modification, and redistribution, including commercial use, subject to its source-availability and other conditions. [IAPWS package license](https://iapws.readthedocs.io/en/latest/modules.html) · [GNU guidance on GPL libraries](https://www.gnu.org/licenses/gpl-faq.en.html#IfLibraryIsGPL)

Third-party components retain their copyright holders and license terms. Preserve the notices below and the license files included with those components when redistributing this project.

## Bundled source code

| Component | Version | Copyright and license | Included records |
| --- | --- | --- | --- |
| [ExerPy](https://github.com/oemof/exerpy) | 0.2.2 | Technische Universität Berlin, Chair of Energy Engineering and Climate Protection, and the ExerPy contributors. MIT. | [License](src/aspen_energy_exergy/_vendor/exerpy/LICENSE), [source manifest](src/aspen_energy_exergy/_vendor/exerpy/vendor_manifest.json), [adapter notice](src/aspen_energy_exergy/_vendor/exerpy/ASPEN_ADAPTER_NOTICE.txt) |
| [TESPy](https://github.com/oemof/tespy) | 0.12.0.post1 | Francesco Witte and the contributors identified in the source headers. MIT. | [License](src/aspen_energy_exergy/_vendor/tespy/LICENSE), [source manifest](src/aspen_energy_exergy/_vendor/tespy/VENDOR_MANIFEST.json) |

The bundled engines include changes that defer optional property-library imports. Their manifests and patch files record the upstream versions and the changes. The original copyright and MIT permission notices remain with the source. These engines are used for balance evaluation with supplied Aspen states. Their numerical equations and thermodynamic methods should be credited to their respective upstream projects.

## Installed dependencies

The package installer obtains the following dependencies separately. The versions below were inspected for this publication. The permitted installation ranges appear in `pyproject.toml`, and each installed distribution carries its applicable license texts and additional bundled-code notices.

| Dependency | Inspected version | Main license | Upstream license or project |
| --- | --- | --- | --- |
| thermo | 0.6.1 | MIT | [Caleb Bell / thermo](https://github.com/CalebBell/thermo/blob/master/LICENSE.txt) |
| chemicals | 1.5.2 | MIT | [Caleb Bell / chemicals](https://github.com/CalebBell/chemicals/blob/master/LICENSE.txt) |
| fluids | 1.3.1 | MIT | [Caleb Bell / fluids](https://github.com/CalebBell/fluids/blob/master/LICENSE.txt) |
| iapws | 1.5.5 | GPLv3 | [Juan José Gómez Romera / iapws](https://github.com/jjgomera/iapws/blob/master/LICENSE) |
| NumPy | 2.5.4 | BSD-3-Clause, with additional bundled-code licenses | [NumPy license](https://github.com/numpy/numpy/blob/main/LICENSE.txt) |
| pandas | 3.1.1 | BSD-3-Clause | [pandas license](https://github.com/pandas-dev/pandas/blob/main/LICENSE) |
| SciPy | 1.18.1 | BSD-3-Clause, with additional bundled-code licenses | [SciPy license](https://github.com/scipy/scipy/blob/main/LICENSE.txt) |
| Pint | 0.26.2 | BSD-3-Clause | [Pint license](https://github.com/hgrecco/pint/blob/master/LICENSE) |
| tabulate | 0.10.0 | MIT | [tabulate license](https://github.com/astanin/python-tabulate/blob/master/LICENSE) |
| msgspec | 0.22.0 | BSD-3-Clause | [msgspec license](https://github.com/msgspec/msgspec/blob/main/LICENSE) |

For redistribution of a full environment, retain the complete license and attribution files supplied by every included distribution, including transitive dependencies and bundled data. For example, thermo includes ChemSep interaction-parameter data with an Artistic License 2.0 notice, and binary NumPy and SciPy distributions include notices for their bundled numerical libraries.

## Scientific methods and references

The elemental exergy reference follows Tang et al. (2021), equation 5 and table 1. Water-property estimates use the Python `iapws` implementation of the IAPWS-IF97 formulation. Mixture-property calculations use thermo, chemicals, and fluids from the Chemical Engineering Design Library. The README describes the assumptions, equations, and limits of these calculations. [Tang et al., DOI 10.3303/CET2183089](https://doi.org/10.3303/CET2183089) · [IAPWS industrial formulation](https://iapws.org/technical-guidance/release/IF97-Rev)

Aspen Plus is a product of Aspen Technology. The project reads user-provided simulation results. Aspen Technology's product names and software remain subject to their respective rights and terms.
