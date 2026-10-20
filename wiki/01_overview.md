# Computational overview

The computational work of BruSafe has two arms. They answer different questions and are connected
at two points, described below.

| Arm | Question | Method | Location |
|---|---|---|---|
| **Molecular dynamics** | Does the single-chain (scCP) linker perturb the coat-protein dimer, and does it change the binding of the dimer to the *pac* hairpin? | all-atom MD with GROMACS (CHARMM36m), followed by MM/PBSA for ΔΔG | `scripts/`, `mdp/`, `analysis/`, `mmpbsa/`, `docs/` |
| **CAPSID framework** | How much correctly formed, cargo-loaded VLP does a culture yield, and which step limits it? | five chained kinetic and statistical models, each calibrated on its own observable | `capsid/` |

## The CAPSID chain

![Overview of the five modules](../capsid/figures/fig28_wiki_model_overview.png)

*Figure 1. The five modules of CAPSID. Each module is fitted to its own observable (bottom row), and
the arrow between two boxes is the only quantity that the next module receives.*

| Module | Computes | Measured in the laboratory by | Notebook |
|---|---|---|---|
| M1 · Expression | coat-protein (CP) concentration over time | biomass (g/L) and CP (mg/L) | 01, 02 |
| M2 · Soluble partition | soluble fraction φ of the CP | soluble and pellet CP (mg/L) | 03 |
| M3 · Capsid assembly | capsid fraction and nucleation lag | SEC-MALS and TEM diameter | 04 |
| M4 · Cargo packaging | cargo copies per particle (from *pac* binding) | dPCR copies and protected fraction | not yet implemented (exploratory notebooks only) |
| M5 · Quantification | VLP concentration at harvest | particles per mL and mg/L of VLP | 05 |

## Connections between the two arms

1. **MM/PBSA → M4.** The binding of the CP dimer to the *pac* hairpin, expressed as a ΔΔG between the
   wild-type and the single-chain dimer (S9 against S4), constrains the effect of the linker on the
   dissociation constant K<sub>d</sub>(pac) that M4 will use.
2. **MD → M3.** The conformer populations of the FG loop from the A/B and C/C dimer simulations
   inform the quasi-equivalence argument behind the *T* = 3 contact schedule of M3.

## Status

| Component | Status |
|---|---|
| CAPSID M1, M2, M3, M5 | implemented; calibrated on synthetic data and confronted with published measurements (pages 03–06) |
| CAPSID M4 | specified; exploratory notebooks only |
| MD campaign | protocol, scripts, structure builds and analysis pipeline complete (page 02); the apo systems S1 and S6 have finished, and the remaining systems are run on the laboratory machines; production results are not part of this repository |

The limitations of every component are collected on [page 09](09_limitations.md).
