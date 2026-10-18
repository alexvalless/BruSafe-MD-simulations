# Limitations and what has not been validated

This page lists what the computational work does not show. Each item names the observation that
would resolve it.

## CAPSID framework

| Limitation | Consequence | What would resolve it |
|---|---|---|
| The contact schedule and degeneracy vector of M3 are placeholders with the correct total (§2 of notebook 04) | the kinetic capsid fraction at 18 h varies between 0.003 and 0.894 across the sensitivity scan; the nucleation lag, the stalled fraction and the intermediate size distribution cannot be quoted | enumerating the contacts from PDB 2MS2 and 1AQ3, which needs no bench work |
| Held-out dose coverage is 50% (IPTG) and 60% (lactose) against a nominal 94% | the model is overconfident at an unseen dose | a Kennedy–O'Hagan discrepancy term, before any forecast interval at an unseen condition is quoted |
| k<sub>deg</sub> has no MS2-CP-specific literature value, and a 5 h time course cannot constrain it | the predicted optimal IPTG dose changes direction with the horizon depending on k<sub>deg</sub> | time points at 12 and 24 h, or a translation-halted arm |
| The leak ε and the gain G are correlated (r ≈ −0.98) in the lactose arm | yield forecasts inherit the ridge | an uninduced arm carried for the whole run |
| The recovery efficiencies of M5 are assumed and not measured | they carry 55% of the variance in cargo copies per particle | the three recovery controls (release, extraction, carryover) |
| M5 absolute figures (harvest titre, doses per litre, cost per dose) use placeholder inputs and an assumed dose size | only the error budget and the factor ranking are quantitative | measured recoveries, and the dose size from the immunology work |
| M4 (CP–*pac* energetics and packaging calibration) has no production notebook | the link between MD and the packaging ratio has not been implemented | the M4 notebook, fed with K<sub>d</sub>(pac) from MM/PBSA |
| The cascade has no off-pathway or malformed species | the yield predictions are optimistic if TEM shows malformed particles that increase with induction | TEM counts across the induction titration |
| The effect of pac binding on nucleation is a qualitative coupling | no quantitative prediction for 0× versus 3× *pac* in vivo | in vitro reassembly at 20–50 µM dimer ± TR hairpin |
| The temperature axis of M2 depends on stage-1 parameters that also change with temperature | one M1 posterior cannot be reused across temperatures | a separate M1 fit for each induction temperature |
| Calibration is on synthetic data; the literature datasets are small (two points for M2, a design without a temperature axis for M1) | the published confrontations test qualitative claims only | the planned wet-lab experiments |

## Molecular dynamics

| Limitation | Consequence |
|---|---|
| Production results for the tiered systems are not part of this repository | the MD pages describe the protocol and the analysis; claims about the linker or about binding are made only after replica spread and convergence have been shown |
| Absolute MM/PBSA binding free energies are unreliable for charged protein–RNA interfaces (typical error 5–20 kcal mol⁻¹) | only ΔΔG between systems is reported, with the sensitivity to the interior dielectric constant |
| Mg²⁺ ions do not equilibrate on the 250 ns timescale | the placement is documented, and the interface check is made before MM/PBSA |
| The force-field choice (CHARMM36m) involves a less standard PB radii setup (`PBRadii = 7`) | if the calibration regression performs poorly, the force field is one of the suspects (`docs/DECISIONS.md`) |
| Review of the force-field decision by the external MD advisor was pending when the decision was logged | to be confirmed in writing |

## Documentation that is not part of the repository

The framework document to which the section numbers of the notebooks refer, the theory compendium
of the exploratory M4 notebooks (`M4_theory_compendium.pdf`), and the bench-data template
(`M4_bench_data_template.csv`) are not included.
