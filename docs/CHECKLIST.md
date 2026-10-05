# BruSafe MD Campaign — Operational Checklist

Model deadline **13 Oct** · Wiki freeze **21 Oct** · Start **7 Sep**

Rule: every system below has a named consumer. If a system has no consumer, cut it.

---

## PHASE 0 — Locks (finish by 9 Sep, then never change)

- [ ] Force field family chosen and written down
      - AMBER ff99SB-ILDN (or ff14SB) + OL3 for RNA → cleanest path for MM/PBSA
      - CHARMM36m → better loop dynamics, slightly less standard PB radii
- [ ] Water model + ion set (TIP3P, 150 mM KCl, Mg²⁺ for RNA-containing systems)
- [ ] GROMACS version pinned, identical on all 10 machines (`gmx --version` logged)
- [ ] Naming convention: `<system>_<conformer>_rep<N>`
- [ ] Storage + backup plan (~40 GB total; offsite copy of every `.xtc` + `.tpr` + `.mdp`)
- [ ] Git repo created; all `.mdp`, `.ndx`, topologies, scripts committed from day one
- [ ] Force-field decision reviewed with Diego — get it in writing for the wiki

## PHASE 1 — Structure prep (by 11 Sep)

- [ ] Extract **C/C symmetric** dimer from 1MSC (not the default chain grab)
- [ ] Extract **A/B asymmetric** dimer separately
- [ ] Strip crystallographic waters, alt-locs, cryoprotectants
- [ ] Protonation states at pH 7.4 (check His tautomers near the RNA face)
- [ ] Build scCP dimer model (AF3 or linker graft onto 1MSC)
      - [ ] Superpose onto 1MSC dimer, report backbone RMSD **before** simulating
      - [ ] Confirm subunit relative orientation is native, not an AF3 artefact
- [ ] Obtain CP–operator cocrystal for all RNA-bound systems (do **not** dock de novo)
- [ ] Build designed pac hairpin by in-place mutation of the cocrystal RNA
- [ ] Build scrambled/non-cognate hairpin control
- [ ] Build literature calibration variants (target 6–8, must include the tighter-binding U(−5)C)
- [ ] Build free hairpin (RNA only, no protein) for the pre-organisation penalty

## PHASE 2 — Equilibration (identical protocol, every system)

- [ ] Rhombic dodecahedron, 1.2 nm padding, solvate, neutralise, 150 mM salt
- [ ] Energy minimisation (steepest descent, Fmax < 1000 kJ/mol/nm)
- [ ] NVT 200 ps, position restraints, 310 K — **separate velocity seed per replica**
- [ ] NPT 1 ns, position restraints
- [ ] NPT 5 ns, unrestrained
- [ ] Gate: density ≈ 1000 kg/m³, T and P stable, no restraint-release blowup

## PHASE 3 — Production

Benchmark first (`bench_4070.sh`), then lock one mdrun config for the whole campaign.

### Tier 1 — cannot ship without these

| # | System | Replicas × length | Consumer |
|---|---|---|---|
| S1 | WT CP dimer, C/C apo | 3 × 250 ns | Reference baseline; all RMSF comparisons |
| S2 | scCP dimer + linker, apo | 3 × 250 ns | **The construct question** — M0/M1 |
| S3 | CP dimer + WT operator | 3 × 250 ns | M3.5 reference complex |
| S4 | CP dimer + designed pac | 3 × 250 ns | M3.5 → K_d prior for M4 |
| S5 | CP dimer + scrambled hairpin | 3 × 250 ns | M3.5 specificity ΔΔG vs dPCR 0×pac arm |

### Tier 2 — do if Tier 1 is clean by 27 Sep

| # | System | Replicas × length | Consumer |
|---|---|---|---|
| S6 | WT CP dimer, A/B apo | 3 × 250 ns | FG-loop switch → M3 quasi-equivalence |
| S7 | Calibration variants (6–8) | 1–2 × 150 ns each | M3.5 method validation regression |
| S8 | Free designed pac hairpin | 3 × 250 ns | Pre-organisation penalty; validates ViennaRNA fold |
| S9 | scCP dimer + designed pac | 3 × 250 ns | Does the linker perturb RNA binding? |

### Tier 3 — stretch only

- [ ] S10: 5-dimer capsid patch (~150k atoms) — scCP in capsid context

### Hard rule
- [ ] **No new systems after 27 Sep.** Anything not launched by then is out of scope.

## PHASE 4 — Analysis (per system)

- [ ] PBC cleanup: `-pbc whole` → `-pbc nojump` → `-center -pbc mol -ur compact`
- [ ] Backbone RMSD (core β-sheet separately from whole chain)
- [ ] RMSF per residue, all replicas overlaid
- [ ] Radius of gyration
- [ ] `gmx dssp` secondary structure timeline
- [ ] Dimer interface: buried SASA, interchain contacts, H-bond occupancy
- [ ] FG-loop specific: loop RMSD, clustering, conformer populations
- [ ] PCA / essential dynamics; **eigenvector overlap S1 vs S2** — the rigorous
      way to claim the linker does not alter the dynamics
- [ ] Convergence evidence:
      - [ ] Block averaging
      - [ ] Cosine content of PC1 (near 1 ⇒ random diffusion, not real motion)
      - [ ] Replica agreement — this is your error bar, not frame-count SEM
- [ ] RNA systems: RNA RMSD, base-pair stability, key-contact H-bond occupancy,
      water-mediated bridges

## PHASE 5 — M3.5 MM/PBSA (start 28 Sep)

- [ ] Strip water/ions from trajectories
- [ ] Per-residue decomposition **first** — must recover known operator contacts.
      Gate: if hot spots miss the crystallographic contacts, stop and fix setup.
- [ ] PB primary, GB cross-check
- [ ] `indi` sensitivity sweep: 1, 2, 4 — report spread
- [ ] Frame-window sensitivity: two disjoint halves
- [ ] Entropy (IE and C2) reported separately, never folded silently into "ΔG"
- [ ] Calibration regression: computed ΔΔG vs published ΔΔG; report Spearman ρ and slope
- [ ] Convert specificity ratios via ΔΔG = −RT ln(ratio), RT = 0.616 kcal/mol at 310 K
- [ ] Hand K_d(pac) with honest uncertainty to M4 as a Bayesian prior

## PHASE 6 — Integration (by 11 Oct)

- [ ] M3.5 → M4: K_d prior wired into the SMC calibration
- [ ] MD → M3: FG-loop conformer populations feed the T=3 assembly argument
- [ ] Prediction registered **before** dPCR readout: expected 3×/1×/0× pac ordering
      and rough magnitude
- [ ] Figure set drafted
- [ ] Explicit "what we did not validate" section written

## PHASE 7 — Reproducibility deliverables

- [ ] Repo public: all mdp, topologies, index files, analysis scripts
- [ ] Systems table: atom counts, box, lengths, replicas, achieved ns/day
- [ ] Hardware + software versions documented
- [ ] Every claim on the wiki traceable to a specific run directory

---

## Do NOT attempt
- Full T=3 capsid all-atom MD (millions of atoms — not on 4070s, not in 5 weeks)
- Coarse-grained assembly MD unless someone already knows Martini
- FEP / umbrella sampling for K_d — no time to do it properly
- Simulating the full kilobase cargo mRNA
- Replica exchange
