# Decision log

This log is append-only. Each entry carries a date and the reason for the decision, so that the
choices made during the campaign (the force field, for example) can be traced and are not reopened
without new evidence.

---

## 2026-09-07 — Force field: CHARMM36m

**Decision.** CHARMM36m parameters for protein and nucleic acids, CHARMM-modified TIP3P water,
150 mM KCl, and Mg²⁺ in RNA-containing systems only.

**Alternative considered.** AMBER ff99SB-ILDN or ff14SB with OL3 for RNA. That combination offers a
more direct path into gmx_MMPBSA: the topologies are native AMBER, the mbondi2 radii are the
validated PB default, and most published protein–RNA MM/PBSA benchmarks use it. CHARMM36m was
chosen in spite of a less standard PB radii setup (`PBRadii = 7`, charmm_radii) because of its
description of loop dynamics, which matter for the FG loop, and because it provides one
self-consistent parameter family for protein, RNA and ions.

**Consequence to track.** `PBRadii = 7` is used in every gmx_MMPBSA run. If the calibration
regression of M3.5 performs poorly, the force field is one of the possible causes, and this entry
is where that should be recorded.

**Review.** Pending: Diego Velasco-González (external MD advisor). Written confirmation was to be
obtained before 2026-09-09, since the choice is a physical judgement and the attribution
strengthens the wiki.

---

## 2026-09-07 — Mg²⁺ handling in RNA systems

**Problem.** Water exchange at Mg²⁺ takes microseconds. On a 250 ns trajectory the ions do not
equilibrate and remain close to where they were placed, so random placement with `genion` adds
noise and no realism.

**Decision.** Placement methods are used in the following order of preference:

1. crystallographic Mg²⁺ sites from the cocrystal, placed explicitly;
2. CHARMM-GUI placement with its ion-distribution options;
3. `genion`, only if 1 and 2 are unavailable, with the placement documented here and the
   Mg²⁺ interface check of `04_postprocess.sh` (`mindist_mg.xvg`) inspected.

**Consequence to track.** An Mg²⁺ ion located in the protein–RNA interface contaminates every M3.5
energy computed from that replica. The check has to be made before MM/PBSA is run.

---

## 2026-09-07 — Replica independence

**Decision.** Velocity seeds are generated for each replica at the NVT stage from a deterministic
hash of `<system>_rep<N>` (see `replica_seed` in `scripts/lib.sh`). Replicas are never branched
from a shared equilibrated frame.

**Reason.** Three replicas that share an equilibration are correlated copies of a single run, and
error bars computed from them would be too small.

---

## 2026-10-03 — Windows (Git Bash) as a supported platform

**Decision.** The pipeline also runs on laboratory machines where WSL cannot be enabled, using
Git Bash and one portable CUDA build of GROMACS copied to every machine. Setup is described in
`docs/WINDOWS.md`.

**Consequence to track.** On Windows, `-pin on` is omitted from the mdrun flags because GROMACS does
not support thread pinning there. Pinning affects speed and not the trajectory physics, so replicas
remain comparable as long as the GROMACS version is identical on every machine (`00_env.sh` records
it). MM/PBSA remains on Linux, WSL or macOS.

---

## Template

## YYYY-MM-DD — <decision>

**Decision.**

**Alternative considered.**

**Consequence to track.**
