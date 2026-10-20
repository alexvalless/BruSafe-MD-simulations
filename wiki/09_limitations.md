# Limitations and what has not been validated

This page lists what the computational work does not show. Each item names the observation that
would resolve it.

## CAPSID framework

| Limitation | Consequence | What would resolve it |
|---|---|---|
| Notebook 04 still runs with a placeholder contact schedule and degeneracy vector that have the correct total (§2 of the notebook); the schedule enumerated from PDB 2MS2 (`capsid/data/structure/`) has not been wired in | the kinetic capsid fraction at 18 h varies between 0.003 and 0.894 across the sensitivity scan; the nucleation lag, the stalled fraction and the intermediate size distribution cannot be quoted | replacing the placeholder with the enumerated schedule and rerunning the notebook. The enumeration is geometric (4.5 Å, at least 5 residue pairs per interface, all interfaces equivalent) and was not repeated on 1AQ3 |
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
| Production results are not part of this repository, and only S1 and S6 had finished when the protocol was written | the MD pages describe the protocol and the analysis; claims about the linker or about binding are made only after replica spread and convergence have been shown |
| Absolute MM/PBSA binding free energies are unreliable for charged protein–RNA interfaces (typical error 5–20 kcal mol⁻¹) | only ΔΔG between systems is reported, with the sensitivity to the interior dielectric constant |
| Mg²⁺ ions do not equilibrate on trajectories of 30–100 ns, and 1ZDH has no crystallographic Mg²⁺ | ions are placed independently for each replica, and the interface table is inspected before MM/PBSA |
| The crystal RNA of the complexes has a mean occupancy of 0.56, only 13 of 19 nucleotides are resolved, and the terminal pair of the stem may fray within 30 ns | the base pairing is checked in every replica; if the stem opens in most of them, a restraint or an extended stem is applied to all systems alike |
| The scCP graft keeps the FG-loop conformation of the second copy, which differs from the RNA-bound chain of S4 (3.07 Å all-Cα RMSD against 0.54 Å for the rigid core) | reported as a limitation of the model in the S9 results |
| The scrambled-RNA control and the calibration variants were dropped, so nothing calibrates the absolute MM/PBSA energies | only ΔΔG for S9 against S4 and for S3 against S4 is reported, and no specificity claim is made |
| Runs of 30 ns per replica sample only the bound basin | the replicas are independent, and the uncertainty is taken from their spread |
| The force-field choice (CHARMM36m) involves a less standard PB radii setup (`PBRadii = 7`) | the force field is one of the possible causes of any anomaly in the MM/PBSA results (`docs/DECISIONS.md`) |
| Review of the force-field decision by the external MD advisor was pending when the decision was logged | to be confirmed in writing |

## Further documentation

The framework document to which the section numbers of the notebooks refer, the results write-up and
the theory compendium of the exploratory M4 notebooks are in `capsid/docs/`.
