# Pre-registered predictions

Each prediction is entered here and committed before the corresponding wet-lab readout exists, so
that the `git log` timestamp documents the order of events.

A pre-registered prediction that agrees with the measurement is stronger evidence than a post-hoc
fit. One that disagrees, reported together with a diagnosis, is still more informative than a model
adjusted after the data have been seen.

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

This is the most informative test. If the pac sites bind independently and packaging requires at
least one nucleation event, the ratio is not 3: near-saturating occupancy drives it toward 1, and
occupancy far from saturation drives it toward 3. The measured ratio therefore locates K_d relative
to the intracellular CP concentration, and M3.5 together with M4 can predict it without any
absolute binding energy being correct.

**Derived from:** M3.5 K_d(pac) = ______ ± ______ , M4 run ______

---

## P2 — Specificity ΔΔG consistency (withdrawn)

**Status:** withdrawn on 2026-10-04. The scrambled-RNA control (S5) and the calibration variants were removed
from the campaign (`docs/DECISIONS.md`), so there is no specificity ΔΔG to compare with the model. The
entry is kept for the record.

Model specificity ratio (P1/P2): ______
Implied ΔΔG = −RT ln(ratio), RT = 0.616 kcal/mol at 310 K: ______ kcal/mol
MM/PBSA ΔΔG(cognate pac − scrambled): ______ ± ______ kcal/mol

**Agreement expected within:** ______ kcal/mol
**If MM/PBSA over-separates** (a known failure mode for charged protein–RNA interfaces, more
pronounced with GB than with PB), this is reported as such and the model is not adjusted to
match.

---

## P3 — Does the linker perturb the dimer?

Predicted essential-subspace overlap, S1_wt_cc_apo vs S2_sccp_apo, first 10 eigenvectors: ______

Predicted interface buried-SASA difference: ______ nm²
Inter-replica spread on that quantity: ______ nm²

The difference counts only if it exceeds the inter-replica spread. Both numbers are written down
before the analysis is inspected.
