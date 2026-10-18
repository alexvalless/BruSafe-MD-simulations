# Molecular dynamics: methods

## Scope

The all-atom campaign addresses two questions about the MS2 coat protein (CP):

1. Does the single-chain linker (scCP) alter the dynamics of the CP dimer relative to the wild-type
   dimer? The comparison is made between **S1** and **S2**.
2. How strongly does the CP dimer bind a designed *pac* hairpin, relative to the wild-type operator
   and to a scrambled hairpin? The comparison uses MM/PBSA and is reported as ΔΔG only.

Simulating the full *T* = 3 capsid, coarse-grained assembly, free-energy perturbation or umbrella
sampling for K<sub>d</sub>, the complete kilobase cargo mRNA, and replica exchange are outside the
scope of the campaign.

## Simulation parameters

| | |
|---|---|
| Force field | CHARMM36m (protein and nucleic acids) |
| Water | CHARMM-modified TIP3P |
| Salt | 150 mM KCl |
| Divalent ions | Mg²⁺ in RNA-containing systems only |
| Temperature | 310 K, V-rescale |
| Pressure | 1 bar, C-rescale, isotropic |
| Time step | 2 fs, constraints on bonds involving hydrogen |
| Box | rhombic dodecahedron, 1.2 nm padding |
| Production | 250 ns per replica (150 ns for the calibration variants), 10 ps frame interval |
| Hardware | 10 workstations with one RTX 4070 each; one fixed GROMACS version and one fixed mdrun configuration |

The choice of CHARMM36m and its consequences (in particular `PBRadii = 7` in gmx_MMPBSA) are
documented in `docs/DECISIONS.md`.

## Systems

| System | Tier | Source | RNA | Replicas × ns | Purpose |
|---|---|---|---|---|---|
| S1 `wt_cc_apo` | 1 | 1MSC, chains C/C | no | 3 × 250 | reference baseline for RMSF and PCA comparisons |
| S2 `sccp_apo` | 1 | model | no | 3 × 250 | does the linker perturb the dimer fold? |
| S3 `cp_operator` | 1 | cocrystal | yes | 3 × 250 | M3.5 reference complex with a published K<sub>d</sub> anchor |
| S4 `cp_pacdesign` | 1 | cocrystal | yes | 3 × 250 | K<sub>d</sub>(pac) prior for M4 |
| S5 `cp_scramble` | 1 | cocrystal | yes | 3 × 250 | specificity ΔΔG against the dPCR 0×pac arm |
| S6 `wt_ab_apo` | 2 | 1MSC, chains A/B | no | 3 × 250 | FG-loop switch, quasi-equivalence in M3 |
| S7 `cal_var01–05`, `cal_u5c` | 2 | cocrystal | yes | 2 × 150 each | MM/PBSA validation regression; `u5c` is a tighter-than-wild-type binder |
| S8 `pac_free` | 2 | model | RNA only | 3 × 250 | pre-organisation penalty; validates the ViennaRNA fold |
| S9 `sccp_pac` | 2 | model | yes | 3 × 250 | does the linker perturb RNA binding? |
| S10 `patch5` | 3 | 2MS2 | no | 1 × 100 | scCP in the capsid context (interface sanity check) |

The registry of systems, with replica numbers and consumers, is `config/systems.tsv`. Every system
has a named consumer, that is, the module or claim it supports.

## Protocol

1. **Preparation.** Solvation, neutralisation and 150 mM salt (`01_prepare.sh`). RNA-containing
   systems are built with CHARMM-GUI Solution Builder, because the `pdb2gmx` residue mapping for RNA
   under the CHARMM36 port is fragile and CHARMM-GUI handles Mg²⁺ placement. The CP–operator
   cocrystal is the starting point for every RNA-bound system; the designed *pac* and the
   calibration variants are obtained by mutating the cocrystal RNA in place.
2. **Equilibration** (`02_equilibrate.sh`): steepest-descent minimisation, NVT for 200 ps and NPT for
   1 ns with position restraints, then NPT for 5 ns without restraints. The velocity seed is drawn
   independently for each replica from a hash of the system and replica name, and replicas are never
   branched from a common equilibrated frame.
3. **Equilibration gate.** Density close to 1000 kg m⁻³, temperature within 2 K of 310 K, and no
   remaining drift of the density between the two halves of the unrestrained run.
4. **Production** (`03_production.sh`): restartable from checkpoints, with a provenance file written
   beside each trajectory (host, GROMACS version, flags, seed, start and end times).
5. **Post-processing** (`04_postprocess.sh`): periodic-boundary correction (whole → no-jump →
   centred compact), backbone RMSD, Rg, per-residue RMSF, SASA, DSSP, and PCA of the Cα atoms. For
   RNA systems: RNA RMSD, hydrogen bonds, and the minimum distance between Mg²⁺ and the protein–RNA
   interface. The first 50 ns of each trajectory are excluded from RMSF and PCA.
6. **MM/PBSA** (M3.5, `run_mmpbsa.sh`): see below.

## Convergence and statistical treatment

- **Error bars come from the spread among independent replicas.** The standard error over correlated
  frames underestimates the uncertainty by roughly an order of magnitude; the factor obtained from
  the campaign's own data (`convergence.py block`) is quoted in the results.
- **Cosine content of PC1** is computed for every replica. A value close to 1 indicates that the
  projection is indistinguishable from random diffusion, and the motion is then not interpreted.
- **Replica agreement** of the per-residue RMSF is evaluated within each system. A difference
  between two systems is considered only if it exceeds the spread among the replicas of each system.
- **S1 versus S2** is compared through the overlap of the essential subspaces over the first 10
  eigenvectors. An overlap above 0.7 supports the statement that the linker does not alter the
  dynamics of the dimer; below 0.5 indicates a real change. The covariance analysis is restricted to
  the CP backbone atoms present in both systems.

## MM/PBSA

- The per-residue decomposition is examined **first**. Its hot spots must recover the known
  operator-contact residues of MS2 CP; otherwise the setup is corrected before any energy is
  interpreted.
- PB is the primary method and GB a consistency check.
- The interior dielectric constant (`indi`) is varied over 1, 2 and 4, and two disjoint frame
  windows are analysed; the resulting spread is part of the error bar.
- Entropy (interaction entropy and C2) is reported separately and never folded into a ΔG.
- Only **ΔΔG between systems** is reported. Absolute MM/PBSA binding free energies are commonly off
  by 5–20 kcal mol⁻¹, and more so for charged protein–RNA interfaces. ΔΔG = −RT ln(ratio), with
  RT = 0.616 kcal mol⁻¹ at 310 K, so that one order of magnitude in K<sub>d</sub> corresponds to
  about 1.4 kcal mol⁻¹.
- The calibration regression (computed against published ΔΔG for the S7 variants) is reported with
  Spearman's ρ and the slope.

## Handling of Mg²⁺

Water exchange at Mg²⁺ takes microseconds, so on a 250 ns trajectory the ions do not equilibrate and
stay close to where they were placed. Crystallographic sites are therefore used where available,
CHARMM-GUI placement otherwise, and `genion` only as a documented last resort. The minimum distance
between each Mg²⁺ ion and the interface is inspected for every RNA-containing replica, because an
ion that sits in the interface contaminates every energy computed from that replica.

## Pre-registered predictions

Predictions for the dPCR readout (the ratio of cargo copies between the 3×pac and 1×pac constructs),
for the consistency between the model specificity and the MM/PBSA ΔΔG, and for the effect of the
linker on the essential subspace are registered in `docs/PREDICTIONS.md` before the corresponding
data exist.

## Reporting pages

`analysis/md_report.py wiki --all` builds one compact interactive page for each system, with all
replicas overlaid and shown as mean ± SD (RMSD, RMSF, Rg, SASA; RNA RMSD and Mg²⁺ distance for RNA
systems) and one trajectory animation. The full set of figures, including DSSP, PCA and block
averaging, is produced with `analysis/md_report.py report <system>`.
