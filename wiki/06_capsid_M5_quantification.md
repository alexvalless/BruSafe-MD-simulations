# CAPSID M5 · Absolute quantification

Notebook: `BruSafe_05_M5_absolute_quantification.ipynb`

## Purpose

The upstream modules produce ratios or closed shells. M5 provides the measurement that converts them
into absolute, dose-relevant quantities: the number of cargo copies per particle,

CPP = (dPCR copies of RNase-protected cargo per mL) / (particles per mL).

CPP gives the absolute packaging fraction (M4 yields only a ratio between constructs), the
concentration of **loaded** VLP (M3 yields closed shells, empty or full), and doses per litre of
culture. It is also the one place at which the dPCR panel and the assembly module measure the same
physical object, so it is a cross-check and not two unrelated assays.

## The numerator: dPCR of RNase-protected cargo

C<sub>meas</sub> = C<sub>true</sub> ε<sub>rel</sub> ε<sub>ext</sub> + B

| Factor | Effect | Control |
|---|---|---|
| ε<sub>rel</sub>, capsid disruption | incomplete release undercounts | known RNA spiked into intact VLPs, disrupted and recovered |
| ε<sub>ext</sub>, extraction recovery | undercounts | external RNA spike through the full protocol |
| B, protection carryover | surviving free RNA overcounts | naked cargo RNA + RNase, no capsid |

B is the only term that biases the result upward, and it is the control most often omitted. dPCR
counting precision is not the limiting factor: with 26 000 partitions the counting CV is below 1%
across the usable occupancy range (0.77% at the optimal occupancy λ = 1.59). With EvaGreen chemistry
cargo and reference cannot be multiplexed by colour, so the between-well variance has to be carried
explicitly when they are run in separate wells.

## The denominator: counting particles

| Method | Suitable at 27 nm? | Comment |
|---|---|---|
| NTA | marginal | practical floor ~30–50 nm for low-index biologicals |
| TRPS | marginal | smallest pores ~40 nm |
| Negative-stain TEM | yes | counting-limited; needs a bead standard |
| SEC-MALS | yes | absolute molar mass; separates loaded from empty |
| **SEC peak + protein assay** | yes | least expensive; CV 7.7% |

Converting the protein mass of the purified capsid peak to particles through the protein mass per
VLP (2.47 MDa) is more precise than counting particles: TEM has a CV of 12.7% for 200 particles and
5.7% for 1000, against 7.7% for the protein assay. TEM is used as an orthogonal check on shape and
malformed fraction, and the protein assay supplies the particle number.

A280 alone cannot be used on a loaded VLP, because RNA absorbs strongly at 280 nm. A two-wavelength
deconvolution (A260/A280) is ill-conditioned in this regime, with a condition number of 43.8: a 1%
absorbance error becomes up to a 44% error in the protein mass. The recommended split is a
RNA-insensitive protein assay (BCA or Bradford) for protein, A260 for RNA, and dPCR for the cargo.
Comparing the A260-derived RNA content with the dPCR cargo count also gives the packaging
specificity, that is, the fraction of packaged RNA that is the construct.

## Error budget

For the placeholder control values used in the notebook, the estimator gives CPP = 0.62 (95%
interval 0.47–0.85, a span of 1.8×). The share of the log-variance is:

| Source | Share |
|---|---|
| extraction recovery | 55.0% |
| particle count (protein assay) | 25.7% |
| capsid disruption | 18.0% |
| protection carryover | 1.0% |
| dPCR counting | 0.3% |

The recovery controls dominate, which is why they have to be measured and not assumed to be 1.
Running them is worth more than any amount of replication without them: the CV falls from 25% with no
controls to 8.8% with three biological replicates and 6.2% with six, and the return flattens beyond
about six replicates.

![dPCR precision](../capsid/figures/fig20_M5_dpcr_precision.png)

*Figure 1. Counting precision of dPCR (CV) against the mean occupancy λ, for 8,500 and 26,000
partitions. The grey band is the usable range of occupancy.*

## Scale-up

The chain from culture volume to doses is

doses per litre = [CP] × f<sub>sol</sub> × f<sub>cap</sub> × f<sub>loaded</sub> × η<sub>purif</sub> / d<sub>mg</sub>.

**Every input in the notebook is a placeholder**, and the dose size in particular has to come from
the immunology and human-practices work. The quantity that does not depend on the placeholders is
the ranking of the factors by their effect on the harvest titre.

![What moves the harvest titre](../capsid/figures/fig21_M5_harvest_titre_tornado.png)

*Figure 2. Share of the variance of log[VLP] at harvest attributable to each factor. Purification
recovery (58%) and the capsid fraction (24%) are the two factors for which additional effort changes
the result.*

With the placeholder inputs, the harvest titre is 15 mg L⁻¹ of loaded VLP protein (95% interval
5.8–30), equivalent to 3.7 × 10¹² particles per mL. The derived planning figures (about 150 doses
per litre, at an assumed 0.10 mg per dose) depend on a dose size that no experiment of ours
constrains, and they are not predictions of potency. Culture cost per dose covers media and
induction only; purification consumables, QC, fill-finish and regulatory costs are not included and
would dominate a real cost of goods.

## Measurement sequence

1. Purify to a clean SEC capsid peak, and check the trace for free dimer and partial shells.
2. Run the three recovery controls.
3. Split the readout: BCA for protein, A260 for RNA, dPCR for cargo.
4. Use at least three biological replicates, and six if the schedule allows.
5. Use TEM as the orthogonal check on shape and malformed fraction.

Cross-checks that must agree: the specificity from this module against the P(3×)/P(0×) enrichment
from M4, and the particle number from protein mass against the particle number from TEM counts. A
disagreement implies a systematic error in one of the assays; a CPP well above 1 copy per particle
implies fragmentation of the cargo or a leaking protection step.
