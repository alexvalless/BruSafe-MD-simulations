# BruSafe MD campaign — operational checklist

Model deadline 13 Oct · wiki freeze 21 Oct · start 7 Sep

Rule: every system listed below has a named consumer, that is, the module or claim it supports.
A system without a consumer is removed.

---

## Phase 0 — Locked choices (to be completed by 9 Sep and not changed afterwards)

- [ ] Force field family chosen and recorded
      - AMBER ff99SB-ILDN (or ff14SB) + OL3 for RNA: the most direct path to MM/PBSA
      - CHARMM36m: better loop dynamics, slightly less standard PB radii
- [ ] Water model and ion set (TIP3P, 150 mM KCl, Mg²⁺ for RNA-containing systems)
- [ ] GROMACS version pinned and identical on all 10 machines (`gmx --version` logged)
- [ ] Naming convention: `<system>_<conformer>_rep<N>`
- [ ] Storage and backup plan (~40 GB in total; offsite copy of every `.xtc`, `.tpr` and `.mdp`)
- [ ] Git repository created; all `.mdp`, `.ndx`, topologies and scripts committed from the start
- [ ] Force-field decision reviewed by the external MD advisor and confirmed in writing for the wiki

## Phase 1 — Structure preparation (by 11 Sep)

- [ ] Extract the **C/C symmetric** dimer from 1MSC (selected deliberately, not by default chain order)
- [ ] Extract the **A/B asymmetric** dimer separately
- [ ] Remove crystallographic waters, alternate locations and cryoprotectants
- [ ] Assign protonation states at pH 7.4 (check His tautomers near the RNA face)
- [ ] Build the scCP dimer model (AF3 or linker graft onto 1MSC)
      - [ ] Superpose onto the 1MSC dimer and report the backbone RMSD **before** simulating
      - [ ] Confirm that the relative orientation of the subunits is native and not an AF3 artefact
- [ ] Obtain the CP–operator cocrystal for all RNA-bound systems (no de novo docking)
- [ ] Build the designed pac hairpin by in-place mutation of the cocrystal RNA
- [ ] Build the scrambled/non-cognate hairpin control
- [ ] Build the literature calibration variants (target 6–8, including the tighter-binding U(−5)C)
- [ ] Build the free hairpin (RNA only, no protein) for the pre-organisation penalty

## Phase 2 — Equilibration (identical protocol for every system)

- [ ] Rhombic dodecahedron, 1.2 nm padding, solvation, neutralisation, 150 mM salt
- [ ] Energy minimisation (steepest descent, Fmax < 1000 kJ/mol/nm)
- [ ] NVT 200 ps with position restraints at 310 K, with a **separate velocity seed per replica**
- [ ] NPT 1 ns with position restraints
- [ ] NPT 5 ns unrestrained
- [ ] Gate: density ≈ 1000 kg/m³, stable T and P, no instability on release of the restraints

## Phase 3 — Production

The GPU benchmark (`bench_gpu.sh`) is run first, and one mdrun configuration is then fixed for the
whole campaign.

### Tier 1 — required

| # | System | Replicas × length | Consumer |
|---|---|---|---|
| S1 | WT CP dimer, C/C apo | 3 × 250 ns | Reference baseline for all RMSF comparisons |
| S2 | scCP dimer + linker, apo | 3 × 250 ns | **The construct question** (M0/M1) |
| S3 | CP dimer + WT operator | 3 × 250 ns | M3.5 reference complex |
| S4 | CP dimer + designed pac | 3 × 250 ns | M3.5 → K_d prior for M4 |
| S5 | CP dimer + scrambled hairpin | 3 × 250 ns | M3.5 specificity ΔΔG against the dPCR 0×pac arm |

### Tier 2 — if Tier 1 is complete by 27 Sep

| # | System | Replicas × length | Consumer |
|---|---|---|---|
| S6 | WT CP dimer, A/B apo | 3 × 250 ns | FG-loop switch → M3 quasi-equivalence |
| S7 | Calibration variants (6–8) | 1–2 × 150 ns each | M3.5 method-validation regression |
| S8 | Free designed pac hairpin | 3 × 250 ns | Pre-organisation penalty; validates the ViennaRNA fold |
| S9 | scCP dimer + designed pac | 3 × 250 ns | Does the linker perturb RNA binding? |

### Tier 3 — optional

- [ ] S10: 5-dimer capsid patch (~150k atoms), scCP in the capsid context

### Cut-off
- [ ] **No new systems after 27 Sep.** Systems not launched by then are out of scope.

## Phase 4 — Analysis (per system)

- [ ] PBC correction: `-pbc whole` → `-pbc nojump` → `-center -pbc mol -ur compact`
- [ ] Backbone RMSD (core β-sheet separately from the whole chain)
- [ ] RMSF per residue, all replicas overlaid
- [ ] Radius of gyration
- [ ] `gmx dssp` secondary-structure timeline
- [ ] Dimer interface: buried SASA, interchain contacts, H-bond occupancy
- [ ] FG loop: loop RMSD, clustering, conformer populations
- [ ] PCA / essential dynamics, with the **eigenvector overlap S1 vs S2**, which is the quantitative
      way to support the statement that the linker does not alter the dynamics
- [ ] Convergence evidence:
      - [ ] Block averaging
      - [ ] Cosine content of PC1 (near 1 indicates random diffusion and not real motion)
      - [ ] Replica agreement, which provides the error bar (not the frame-count SEM)
- [ ] RNA systems: RNA RMSD, base-pair stability, key-contact H-bond occupancy, water-mediated bridges

## Phase 5 — M3.5 MM/PBSA (from 28 Sep)

- [ ] Strip water and ions from the trajectories
- [ ] Per-residue decomposition **first**; it must recover the known operator contacts.
      Gate: if the hot spots miss the crystallographic contacts, the setup is corrected before
      anything else is done.
- [ ] PB as the primary method, GB as a cross-check
- [ ] `indi` sensitivity sweep at 1, 2 and 4, with the spread reported
- [ ] Frame-window sensitivity: two disjoint halves
- [ ] Entropy (IE and C2) reported separately and never folded silently into "ΔG"
- [ ] Calibration regression: computed ΔΔG against published ΔΔG, with Spearman ρ and slope
- [ ] Specificity ratios converted with ΔΔG = −RT ln(ratio), RT = 0.616 kcal/mol at 310 K
- [ ] K_d(pac) passed to M4 as a Bayesian prior, with its uncertainty

## Phase 6 — Integration (by 11 Oct)

- [ ] M3.5 → M4: K_d prior connected to the SMC calibration
- [ ] MD → M3: FG-loop conformer populations feed the T=3 assembly argument
- [ ] Prediction registered **before** the dPCR readout: expected 3×/1×/0× pac ordering and
      approximate magnitude
- [ ] Figure set drafted
- [ ] Section on what was not validated written explicitly

## Phase 7 — Reproducibility deliverables

- [ ] Repository public: all mdp files, topologies, index files and analysis scripts
- [ ] Systems table: atom counts, box, lengths, replicas, achieved ns/day
- [ ] Hardware and software versions documented
- [ ] Every claim on the wiki traceable to a specific run directory

---

## Out of scope

- Full T=3 capsid all-atom MD (millions of atoms; not feasible on RTX 4070 GPUs within five weeks)
- Coarse-grained assembly MD, unless someone on the team already has Martini experience
- FEP or umbrella sampling for K_d (insufficient time to do them properly)
- Simulation of the full kilobase cargo mRNA
- Replica exchange
