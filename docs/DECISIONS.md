# Decision log

Append-only. Every entry gets a date and a reason. This file is what you point
at when a judge asks "why CHARMM36m?" and it is also what stops the team
quietly re-litigating a settled choice in week four.

---

## 2026-09-07 — Force field: CHARMM36m

**Decision.** CHARMM36m protein + nucleic acid parameters, CHARMM-modified
TIP3P water, 150 mM KCl, Mg²⁺ in RNA-containing systems only.

**Alternative rejected.** AMBER ff99SB-ILDN / ff14SB + OL3. That combination
has the smoother path into gmx_MMPBSA (native AMBER, mbondi2 radii are the
validated PB default, and most published protein–RNA MM/PBSA benchmarks use
it). We accept a slightly less standard PB radii setup (`PBRadii = 7`,
charmm_radii) in exchange for CHARMM36m's loop dynamics, which matter for the
FG loop, and for one self-consistent parameter family across protein, RNA and
ions.

**Consequence to track.** `PBRadii = 7` in every gmx_MMPBSA run. If the
calibration regression in M3.5 comes out poor, the force field is one of the
suspects and this entry is where that gets recorded.

**Reviewed by.** _pending — Diego Velasco-González (external MD advisor)._
Get this in writing before 2026-09-09; it is a physics judgement call and the
attribution strengthens the wiki.

---

## 2026-09-07 — Mg²⁺ handling in RNA systems

**Problem.** Mg²⁺ water-exchange is a microsecond process. On a 250 ns
trajectory the ions do not equilibrate — they stay approximately wherever they
were placed. Random `genion` placement therefore adds noise, not realism.

**Decision.** In priority order:
1. crystallographic Mg²⁺ sites from the cocrystal, placed explicitly;
2. CHARMM-GUI placement using its ion-distribution options;
3. `genion` only if 1 and 2 are unavailable, and only with the placement
   documented here and the `mindist_mg` interface check run in
   `04_postprocess.sh`.

**Consequence to track.** Any Mg²⁺ that parks in the protein–RNA interface
contaminates every M3.5 energy from that replica. Check before running
MM/PBSA, not after.

---

## 2026-09-07 — Replica independence

**Decision.** Velocity seeds are generated per replica at the NVT stage, from
a deterministic hash of `<system>_rep<N>` (see `replica_seed` in
`scripts/lib.sh`). Replicas are never branched from a shared equilibrated
frame.

**Reason.** Three replicas sharing an equilibration are one run in a
trenchcoat, and the error bars computed from them are fiction.

---

## 2026-10-04 — Shorter campaign: HMR at 4 fs, 100 ns apo / 5 × 30 ns complexes

**Problem.** Nine days to the model deadline, and production runs on shared
school machines that are not ours around the clock. 3 × 250 ns per system no
longer fits.

**Decision.**
1. *Lengths follow the claim, not a uniform number.* Apo systems that carry
   PCA / RMSF comparisons (S1, S2, S6): 3 × 100 ns. RNA complexes that feed
   MM/PBSA ΔΔG (S3–S5, S9): 5 × 30 ns; calibration variants (S7): 3 × 30 ns.
   For MM/PBSA, many short independent replicas give better-converged and more
   honest error bars than one long trajectory (ensemble approach, e.g. ESMACS,
   Wan et al.), because the uncertainty is dominated by replica-to-replica
   spread. PCA overlap needs longer single trajectories, which is why the apo
   systems keep 100 ns.
2. *Hydrogen mass repartitioning, dt = 4 fs.* Hydrogen masses ×3 (3.024 Da),
   taken from the bonded heavy atom; water untouched; h-bond constraints as
   before. Written into the topology by `scripts/hmr_top.py` because
   `mass-repartition-factor` only exists from GROMACS 2024. NVT heating stays at
   2 fs. grompp prints a NOTE that some heavy-atom bonds oscillate faster than
   10 × dt; this is expected with HMR at 4 fs and is a note, not a warning.
3. *Box padding 1.0 nm* (was 1.2). Periodic images stay ≥ 2.0 nm apart, above
   the 1.2 nm cutoff. S8 (free hairpin) keeps 1.2 nm (`BRUSAFE_BOX_PAD=1.2`).
   CHARMM-GUI builds use its octahedral box with 10 Å edge distance.
4. *Unrestrained NPT 2 ns* (was 5). The first `skip_ns` of every production run
   (registry column: 10 ns apo, 5 ns complexes) is discarded before any analysis,
   RMSF, PCA and MM/PBSA alike.

**Alternative rejected.** A uniform 20 ns for every system. Cheap, but PCA
subspace overlap does not converge on 20 ns, and "limited compute" is not a
justification a judge accepts; convergence evidence is.

**Consequence to track.** Every claim reports its convergence check
(`convergence.py block`, first-half vs second-half of the analysed window). Any
replica started under the old 2 fs / 1.2 nm protocol is rerun, never mixed with
the new ones.

---

## 2026-10-04 — Cocrystal: 1ZDH, chains A + B(op5) + R

**Decision.** All CP–RNA systems start from 1ZDH (MS2 capsid soaked with the
operator hairpin, 2.7 Å). `00_build_inputs.py show 1ZDH` gives:

- The deposited A, B, C are an asymmetric unit around the quasi-three-fold,
  not a dimer (A–B: 33 contacts). The A/B dimer is A plus B moved by BIOMT
  op5 (171 contacts); `build` picks that operator automatically.
- RNA chain R sits on that A/B dimer (19 + 3 contacts). Chain S overlaps its
  own two-fold copy on the C/C dimer: a symmetry-averaged hairpin, not used.
- Only residues 4–16 of the 19-mer are resolved (13 nt, four base pairs, the
  A6 bulge and the loop). No crystallographic Mg²⁺ or other HETATM.
- The crystallised RNA is the **C(−5) variant** (loop AUCA), not the wild type.
  S7_cal_u5c therefore uses it as is, and S3 (wild type) reverts residue 11 to U.
- Chain R is refined at mean occupancy 0.56 (incomplete soaking of the capsid
  sites). It is the only copy on an A/B dimer (S is symmetry-averaged), so it is
  used anyway: the pose is the bound one, its coordinates are just less certain
  than the protein's. Restrained NVT/NPT lets it settle, and the MM/PBSA
  decomposition gate (known operator-contact residues must come out as hot
  spots) is the check that the interface survived.

**Alternative rejected.** Docking a hairpin onto 1MSC. A solved complex exists.

**Consequence to track.** The 4-bp stem ends in a U–A pair and may fray within
30 ns; check base pairing of U4–A16 per replica. If it opens in most replicas,
either extend the stem by three ideal A-form pairs (all systems alike) or add a
weak flat-bottom restraint on that pair (all systems alike), and record it here.
Every designed sequence (pac, scramble, calibration) must be written over the
same 13-nt frame.

---

## 2026-10-04 — One RNA for every system; S3, S5, S7 dropped

**Decision.** Every RNA-bound system carries the pac hairpin
ACAUGAGGAUCACCCAUGU (the C(−5) operator). Residues 4–16 are exactly the RNA
resolved in 1ZDH chain R, so S4 uses the crystal RNA with no mutation, and only
those 13 nt are simulated: the three terminal pairs (A1–U19, C2–G18, A3–U17)
are not resolved and do not contact the protein. S8 is the S4 hairpin without
protein (bound pose); S9 grafts the scCP model onto S4 and keeps its RNA.

The wild-type operator (S3), the scrambled control (S5) and the calibration
variants (S7) are removed from the registry; S7_cal_u5c would have been an exact
duplicate of S4. This supersedes the S3/S7 sequences in the 1ZDH entry above.

**Consequence to track.** No specificity ΔΔG (P2 in PREDICTIONS.md) and no
calibration regression against literature K_d: nothing validates absolute
MM/PBSA numbers, so they stay unreported (README rule). The one ΔΔG left is
S9 − S4, the linker effect on binding, which makes S9 the only RNA result that
carries a claim. The per-residue decomposition gate on S4 is the remaining check
that the interface is right.

---

## TEMPLATE

## YYYY-MM-DD — <decision>

**Decision.**

**Alternative rejected.**

**Consequence to track.**
