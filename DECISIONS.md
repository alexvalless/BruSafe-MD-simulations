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

## TEMPLATE

## YYYY-MM-DD — <decision>

**Decision.**

**Alternative rejected.**

**Consequence to track.**
