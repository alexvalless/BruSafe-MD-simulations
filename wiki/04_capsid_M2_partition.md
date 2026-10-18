# CAPSID M2 · Soluble partition

Notebooks: `BruSafe_03_M2_soluble_partition.ipynb`,
`literature/BruSafe_03L_M2_literature_arenin.ipynb`

## The gap that M2 closes

M1 predicts the cumulative coat protein (CP). The assembly module (M3) needs assembly-competent
dimer. At T7 flux in BL21(DE3), not every translated chain folds: nascent chains partition between
folding, aggregation into inclusion bodies and proteolysis, and the split depends on the synthesis
flux. M2 supplies the conversion.

## Model

A nascent chain U can fold (k<sub>f</sub>), be degraded (k<sub>d</sub>), or aggregate with an order n
(k<sub>a</sub>). Because these processes (10⁰–10² s) are much faster than expression (10³–10⁴ s), U
is at quasi-steady state, and the fraction φ of chains that escape aggregation is a closed-form
function of the synthesis flux s. For n = 2,

φ(s) = 2(√(1+σ) − 1)/σ, with σ = 8s/s<sub>50</sub>,

with φ → 1 − s/2s<sub>50</sub> at low flux and φ → √(s<sub>50</sub>/2s) at high flux. The temperature
dependence enters through an Arrhenius closure of s<sub>50</sub> with a single parameter, ΔE.

![Partition function](../capsid/figures/fig15_M2_partition_function.png)

*Figure 1. Left: φ against flux at three temperatures; the dotted line is the high-flux asymptote
√(s₅₀/2s). Right: effect of the aggregation order n on the high-flux slope. n is treated as a
sharpness parameter and not as a measured stoichiometry.*

## Bridge from M1

The partition acts on the **flux** of nascent chains and not on the accumulated pool. The pools of
soluble, aggregated and capsid-incorporated CP sum exactly to the M1 pool P(t), for any φ and any
assembly sink, provided that all three carry the same loss rate (growth dilution plus
k<sub>deg</sub>). The notebook verifies this identity numerically (to ~10⁻⁷ of the peak P). It
means that stage 1 (M1 on total CP) is unaffected by stage 2, so that the two-stage fit is
legitimate. Two further results follow:

- Because the soluble lane holds soluble plus capsid CP, the soluble fraction measured on a gel is
  invariant to the assembly sink (to 10⁻⁹). Parameters of M2 can be fitted without M3 and M3 added
  afterwards. The invariance is lost if VLPs pellet in the clearing spin (a shift of 0.50 when they
  do), so the spin speed has to be at most 20 000 g and recorded.
- Giving inclusion bodies their own degradation rate breaks the identity (by about 13%), so that
  variant is only a sensitivity analysis.

## Observation model

The gel reports two concentrations (soluble and pellet CP, in mg L⁻¹), from which the soluble
fraction is derived. It is modelled with a Beta likelihood whose concentration is estimated.
Observations at the boundary (a gel reading a soluble fraction of exactly 1 means "no pellet band
above the detection limit") are treated as **censored** observations.

## Which experimental axis is informative

Across the four-dose lactose panel, the peak flux varies by only about 2.6×, because LacY uptake is
saturated at every dose, and the soluble fraction at 18 h varies by 0.10–0.12, about twice the
densitometry noise. Temperature moves the same quantity by about 0.35 between 37 and 18 °C. The
experiment should therefore be titrated on temperature and not on dose. Temperature also changes
μ<sub>max</sub>, k<sub>cat</sub> and G, so M1 has to be fitted separately for each induction
temperature.

## Calibration study

Recovery of a single synthetic dataset gives one pass/fail outcome per parameter, which cannot
distinguish a biased posterior from an unlucky draw. The module is therefore checked by
simulation-based calibration: parameters are drawn from the prior, data are simulated, the model is
refitted, and the rank of the true value in the posterior is recorded. The ranks should be uniform.

![Calibration ranks](../capsid/figures/fig27_M2_calibration_ranks.png)

*Figure 2. Rank histograms of the true parameter values in the posterior (simulation-based
calibration of M2, 16 replicates). Uniform ranks indicate a calibrated posterior.*

This check identified a defect in the original observation model: 14% of simulated gels returned a
soluble fraction of exactly 1, and the likelihood evaluated a Beta density at the boundary instead
of treating those observations as censored. Before the correction, coverage was 0.64–0.73 against a
nominal 0.94, and the posteriors were about 1.8× too narrow for every parameter. A conventional
recovery check had not detected it. On synthetic data the corrected model recovers 4 of 4 parameters.

## Held-out temperature

The model is fitted at 37 and 25 °C and predicts 18 °C with no free parameter. The Arrhenius closure
has one parameter, so this is an out-of-sample test of the model structure and not an interpolation.
Coverage is computed from the posterior predictive, including the Beta observation noise.

![Held-out temperature](../capsid/figures/fig18_M2_heldout_temperature.png)

*Figure 3. Prediction at the held-out temperature of 18 °C.*

## Confrontation with published data (notebook 03L)

Hidalgo-Vázquez et al. (2025) measured soluble and insoluble arenin in BL21(DE3) at 37 and 30 °C
(0.5 mM IPTG, harvest 8 h after induction): soluble fractions of 0.341 ± 0.024 and 0.686 ± 0.039,
respectively, and a ratio of 1.47 in total target protein between the two temperatures. With only
two measurements and two free parameters, the posterior is a ridge. Two results are independent of
the unknown absolute flux scale:

- The ratio of the s<sub>50</sub> values (37 °C over 30 °C) is 0.173, which gives
  ΔE/R ≈ 23 600 K from the closed-form inversion (the SMC posterior median is ≈ 20 000 K).
- If the change in flux were the only effect of temperature, the soluble fraction at 37 °C would be
  predicted to be 0.62, which is 11.6 measurement standard deviations above the measured 0.34.
  Temperature therefore acts on the partition itself and not only on the flux.

The calibrated closure predicts soluble fractions of 0.87 at 25 °C and 0.97 at 18 °C. One gel at
either temperature would falsify the Arrhenius closure.

## Experimental recommendations

- Load OD-matched samples, with the whole pellet resuspended to the volume of the cleared lysate.
- Use whole-cell lysate for the total-CP band of stage 1.
- Clear at 20 000 g or less, and record the speed.
- Fit M1 separately for each induction temperature.
- Add a lactose arm below K<sub>Y</sub> (100–200 µM) if dose is to be informative.
