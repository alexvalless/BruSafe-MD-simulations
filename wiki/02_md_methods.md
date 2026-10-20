# Molecular dynamics: methods

## Scope

The all-atom campaign addresses two questions about the MS2 coat protein (CP):

1. Does the single-chain linker (scCP) alter the dynamics of the CP dimer relative to the wild-type
   dimer? The comparison is made between **S1** and **S2**.
2. Does the linker alter the binding of the CP dimer to the *pac* hairpin? The comparison uses MM/PBSA
   and is reported as a ΔΔG between **S9** (scCP + pac) and **S4** (wild-type CP + pac).

A second, smaller comparison (**S3** against **S4**) differs in a single nucleotide of the RNA, U or C at
position −5 of the operator.

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
| Time step | 4 fs, constraints on bonds involving hydrogen, hydrogen mass repartitioning (H → 3.024 Da, written into the topology); equilibration heating at 2 fs |
| Box | rhombic dodecahedron, 1.0 nm padding (1.2 nm for the free hairpin) |
| Production | apo systems 3 × 100 ns; RNA complexes 5 × 30 ns; free hairpin 3 × 50 ns; the first 10 ns (apo) or 5 ns (complexes) of each run are discarded before any analysis |
| Hardware | RTX 4070 workstations of the laboratory, one run per GPU; one fixed GROMACS version (2026.3 with CUDA on the laboratory machines) and one fixed mdrun configuration |

The lengths follow the claim each system supports. Principal-component comparisons need long single
trajectories, which is why the apo systems keep 100 ns. For MM/PBSA, many short independent replicas
give better-converged and more reliable error bars than a single long trajectory, since the uncertainty
is dominated by the spread among replicas, and so the RNA complexes use five replicas of 30 ns.
Hydrogen mass repartitioning is applied by `scripts/hmr_top.py` because the corresponding mdp option
exists only from GROMACS 2024. The choice of CHARMM36m and the other decisions of the campaign, with
their reasons, are documented in `docs/DECISIONS.md`.

## Systems

| System | Tier | Source | RNA | Replicas × ns | Purpose |
|---|---|---|---|---|---|
| S1 `wt_cc_apo` | 1 | 2MS2, chains C/C | no | 3 × 100 | reference baseline for RMSF and PCA comparisons |
| S2 `sccp_apo` | 1 | AlphaFold3 model | no | 3 × 100 | does the linker perturb the dimer fold? |
| S4 `cp_pacdesign` | 1 | 1ZDH, chains A, B (operator 5) and R | yes | 5 × 30 | pac complex; MM/PBSA reference for the linker ΔΔG of S9 |
| S3 `cp_operator` | 2 | S4 with residue 11 reverted to U | yes | 5 × 30 | wild-type operator (U at −5) against S4 (C at −5) |
| S6 `wt_ab_apo` | 2 | 2MS2, chains A/B | no | 3 × 100 | FG-loop switch, quasi-equivalence in M3 |
| S8 `pac_free` | 2 | S4 hairpin, bound pose | RNA only | 3 × 50 | pre-organisation penalty |
| S9 `sccp_pac` | 2 | scCP model grafted on S4 | yes | 5 × 30 | does the linker perturb RNA binding? ΔΔG against S4 |
| S10 `patch5` | 3 | 2MS2 | no | 1 × 50 | scCP in the capsid context (interface sanity check) |

The registry, with replica numbers and consumers, is `config/systems.tsv`; the RNA sequence of each
RNA-bound system is in `config/rna_designs.tsv`. Every system has a named consumer, that is, the module
or claim it supports. A scrambled-RNA control and a set of calibration variants were removed on
2026-10-04, so no specificity ΔΔG and no calibration regression against literature K<sub>d</sub> values
are reported.

## Starting structures

- **Apo dimers (S1, S6).** The asymmetric unit of the 2MS2 capsid has chains A, B and C. The A/B dimer
  is obtained from chains A and B, and the C/C dimer from chain C and its two-fold symmetry mate
  (BIOMT operator 6, 154 contacts). `scripts/00_build_inputs.py` selects each second chain with the
  operator that gives the most contacts with what is already assembled.
- **RNA complexes (S3, S4, S8, S9).** The cocrystal 1ZDH (2.7 Å) provides the protein–RNA complex: the
  A/B dimer is chain A plus chain B moved by operator 5, and the RNA is chain R. Only residues 4–16 of
  the 19-mer are resolved, so the 13 nucleotides UGAGGAUCACCCA are simulated: the pac hairpin of S4 is
  the crystal RNA, which is the C(−5) variant, and S3 reverts residue 11 to U. The occupancy of chain R
  is 0.56 on average (incomplete soaking), and no magnesium ion was resolved.
- **scCP (S2, S9).** One chain of 262 residues, CP(1–129) + GGSG + CP(134–262), predicted with
  AlphaFold3 (pLDDT 88.8; 51–63 in the linker). For S9 the model is superposed on the two CP chains of
  S4, and the RNA and Mg²⁺ of S4 are taken over. The rigid-core RMSD of the graft is 0.54 Å for 225 Cα
  atoms; the all-Cα RMSD of 3.07 Å comes entirely from the FG loop of the second copy.
- **Topologies.** S1–S4, S6, S8 and S9 are built with CHARMM-GUI Solution Builder. Only the CHARMM-GUI
  topology is kept: the box, water and ions are rebuilt with GROMACS, because the octahedral box written
  by CHARMM-GUI contained periodic-image clashes in GROMACS. The step-by-step procedure is in
  `docs/TOPOLOGIES.md`.

## Protocol

1. **Preparation** (`01_prepare.sh`). Solvation, neutralisation and 150 mM KCl, hydrogen mass
   repartitioning, and index groups. Mg²⁺ ions are placed with `genion` using a seed derived from the
   system and replica names, so every replica has its own ion placement.
2. **Equilibration** (`02_equilibrate.sh`): steepest-descent minimisation, NVT for 200 ps and NPT for
   1 ns with position restraints (1000 kJ mol⁻¹ nm⁻²), then NPT for 2 ns without restraints. The
   velocity seed is drawn independently for each replica from a hash of the system and replica name, and
   replicas are never branched from a common equilibrated frame.
3. **Equilibration gate.** Density close to 1000 kg m⁻³, temperature within 2 K of 310 K, and no
   remaining drift of the density between the two halves of the unrestrained run. The drift criterion
   is a change of more than 2 kg m⁻³ that also exceeds three block-averaged standard errors, so that the
   noise of a small box is not mistaken for drift.
4. **Production** (`03_production.sh`, `run_slot.sh`): restartable from checkpoints, with a provenance
   file written beside each trajectory (host, GROMACS version, flags, seed, start and end times).
5. **Post-processing** (`04_postprocess.sh`): the solute is gathered with `-pbc cluster` after making
   molecules whole, and the first clustered frame is the reference for RMSD, RMSF and PCA. Backbone
   RMSD, Rg, per-residue RMSF, SASA, DSSP and PCA of the Cα atoms follow, and, for RNA systems, RNA
   RMSD, protein–RNA hydrogen bonds and the Mg²⁺ interface table. The first `skip_ns` of each
   trajectory are excluded from every statistic.
6. **MM/PBSA** (M3.5, `run_mmpbsa.sh`): see below.

## Convergence and statistical treatment

- **Error bars come from the spread among independent replicas.** The standard error over correlated
  frames underestimates the uncertainty by roughly an order of magnitude; the factor obtained from the
  campaign's own data (`convergence.py block`) is quoted in the results.
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
- Only **ΔΔG between systems** is reported, namely S9 against S4 (effect of the linker) and, as a
  one-nucleotide comparison, S3 against S4. Absolute MM/PBSA binding free energies are commonly off by
  5–20 kcal mol⁻¹, and more so for charged protein–RNA interfaces, and nothing in the campaign
  calibrates them. ΔΔG = −RT ln(ratio), with RT = 0.616 kcal mol⁻¹ at 310 K, so that one order of
  magnitude in K<sub>d</sub> corresponds to about 1.4 kcal mol⁻¹.

## Handling of Mg²⁺

Water exchange at Mg²⁺ takes microseconds, so on trajectories of 30–100 ns the ions do not equilibrate
and stay close to where they were placed. Because 1ZDH has no crystallographic Mg²⁺ site, the ions are
placed with `genion`, independently for each replica so that replicas do not share an ion placement.
For every RNA-containing replica, `mg_interface.txt` reports the fraction of frames in which each ion
lies within 0.5 nm of the RNA, of the protein and of both; an ion that sits in the interface is flagged
in the MM/PBSA table and not averaged silently. The Mg²⁺ concentration of the free hairpin S8 is
matched to that of S4 (one ion in the smaller box).

## Pre-registered predictions

Predictions for the dPCR readout (the ratio of cargo copies between the 3×pac and 1×pac constructs)
and for the effect of the linker on the essential subspace are registered in `docs/PREDICTIONS.md`
before the corresponding data exist. The prediction on specificity was withdrawn together with the
scrambled control.

## Reporting pages

`analysis/md_report.py wiki --all` builds one compact interactive page for each system, with all
replicas overlaid and shown as mean ± SD (RMSD, RMSF, Rg, SASA; RNA RMSD and Mg²⁺ distance for RNA
systems) and one trajectory animation. The full set of figures, including DSSP, PCA and block
averaging, is produced with `analysis/md_report.py report <system>`.
