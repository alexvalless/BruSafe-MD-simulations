# Data and sources

The CAPSID notebooks run on synthetic data with known truth until the wet-lab data are available.
In addition, three published datasets are used to confront the models with measurements. The data
tables are in `capsid/data/`, and a data dictionary is in `capsid/data/README.md`.

| Module | Quantity | Source | Locator | Conversion applied |
|---|---|---|---|---|
| M1 | yield against temperature × time × IPTG, 15-run Box–Behnken design, three proteins, three centre-point replicates | Gutiérrez-González M, et al. *Sci Rep* 2019; doi:10.1038/s41598-019-53200-7 (PMC6856375) | Table 2 (levels resolved from the Table 1 legend) | transcribed; protein-specific factor levels expanded to actual values |
| M1 | fed-batch volumetric production at 22, 25 and 30 °C and 6, 8 and 10 h, with OD₆₀₀ and DCW | Carnicer M, et al. *Sci Rep* 2026; doi:10.1038/s41598-026-54839-9 (PMC13458553) | Tables 4 and 5 | transcribed; reference only, because the fed-batch regime is outside the shake-flask assumptions of the model |
| M2 | soluble and insoluble target protein (arenin) at 30 and 37 °C with replicate SD | Hidalgo-Vázquez E, et al. *Bioresour Bioprocess* 2025; doi:10.1186/s40643-025-00986-2 (PMC12681508) | Table 1 | soluble fraction = soluble/(soluble + insoluble); SD propagated by the delta method |
| M3 / M4 | thermal denaturation midpoints of wild-type and mutant MS2 VLPs | Pistono PE, et al. *Biochemistry* 2024; doi:10.1021/acs.biochem.4c00225 (PMC11308365) | Results, SYPRO Orange thermal shift assay | T<sub>m</sub> inverted to ΔG<sub>contact</sub> with Eq. (4.8) at the stated assay concentration |

## Observables expected from the laboratory

The models use physical quantities and not instrument readings:

| Observable | Definition |
|---|---|
| biomass | gDCW L⁻¹, from OD₆₀₀ × 0.40 or from wet pellet weight × 0.25 |
| coat protein | mg L⁻¹, by densitometry against a purified standard curve, one slope per gel |
| soluble / insoluble CP | two concentrations in mg L⁻¹, and not a band ratio |
| correctly formed fraction | from TEM: of at least 200 particles, those with a diameter of 27 ± 2 nm, with a binomial interval |
| particle concentration | from SEC with a protein assay |
| cargo copies | by dPCR of RNase-protected cargo |

## Key model constants and their origin

| Constant | Value | Origin |
|---|---|---|
| CP monomer mass | 13.7 kDa | Golmohammadi et al. 1993 |
| dimers per *T* = 3 capsid | 90 | 60 A/B + 30 C/C dimers; 180 monomers |
| dimer–dimer contacts in the closed shell | 180 | geometry of the *T* = 3 shell |
| total cytoplasmic protein, in CP-mass equivalents | 18 000 µM (~250 mg mL⁻¹) | Bremer & Dennis 2008; Milo et al. 2010 |
| biomass per OD₆₀₀ unit | 0.40 gDCW L⁻¹ | E. coli in rich medium, 0.3–0.5; to be calibrated for the strain in use |

The complete provenance of each equation term (class, citation) is tabulated in notebooks 01 and 02.
