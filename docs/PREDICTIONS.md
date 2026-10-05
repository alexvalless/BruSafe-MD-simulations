# Pre-registered predictions

Fill this in and **commit it before the corresponding wet-lab readout exists**.
Then `git log` timestamps it.

A pre-registered prediction that comes out right is worth more than ten
post-hoc fits. A pre-registered prediction that comes out wrong, reported
honestly with a diagnosis, is still a stronger Best Model entry than a model
quietly tuned to match data it had already seen.

---

## P1 — dPCR packaging ratios across the pac panel

**Registered:** _YYYY-MM-DD, commit ____________
**Readout not yet available:** ☐ confirmed

Predicted ordering of RNase-protected cargo copies:

| Construct | Predicted, relative to MEV+0×pac | 95% interval |
|---|---|---|
| MEV+3×pac | | |
| MEV+1×pac | | |
| MEV+0×pac | 1.0 (reference) | — |
| eGFP+3×pac | | |

**Predicted 3× / 1× ratio:** ______

This one is the sharp test. If pac sites bind independently and packaging
needs at least one nucleation event, the ratio is *not* 3. Near-saturating
occupancy drives it toward 1; far-from-saturating drives it toward 3. So the
measured ratio localises K_d relative to the intracellular CP concentration —
and that is a prediction M3.5 + M4 can make with no absolute energy being
correct.

**Derived from:** M3.5 K_d(pac) = ______ ± ______ , M4 run ______

---

## P2 — Specificity ΔΔG consistency

Model specificity ratio (P1/P2): ______
Implied ΔΔG = −RT ln(ratio), RT = 0.616 kcal/mol at 310 K: ______ kcal/mol
MM/PBSA ΔΔG(cognate pac − scrambled): ______ ± ______ kcal/mol

**Agreement expected within:** ______ kcal/mol
**If MM/PBSA over-separates** (a known failure mode for charged protein–RNA
interfaces, worse with GB than PB): say so plainly rather than adjusting the
model to match.

---

## P3 — Does the linker perturb the dimer?

Predicted essential-subspace overlap, S1_wt_cc_apo vs S2_sccp_apo,
first 10 eigenvectors: ______

Predicted interface buried-SASA difference: ______ nm²
Inter-replica spread on that quantity: ______ nm²

**The difference only counts if it exceeds the spread.** Write both numbers
down before you look.
