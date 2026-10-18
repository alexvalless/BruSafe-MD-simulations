# Notebook catalogue

All notebooks, grouped by role. The same information in machine-readable form is in
[`../MANIFEST.csv`](../MANIFEST.csv). Notebooks are in English unless the language column says
otherwise.

Each notebook is run from its own folder (see the main [README](../README.md#running-the-notebooks)),
and the files they exchange are written to that folder. Those files (`*.npz`, `*.nc`, and the
posterior CSV files) are excluded from version control.

## 1. Production notebooks

| # | Notebook | Module | Reads | Writes |
|---|---|---|---|---|
| 01 | `BruSafe_01_M1_expression_IPTG.ipynb` | M1 | `VA_timecourse.csv` if `USE_SYNTHETIC = False` | `M1_iptg_posterior_draws.csv`, `M1_iptg_posterior_moments.csv`, `M1_iptg_idata.nc` |
| 02 | `BruSafe_02_M1_expression_lactose.ipynb` | M1 | `VA_lactose_timecourse.csv` if `USE_SYNTHETIC = False` | `M1_lactose_posterior_draws.csv`, `M1_lactose_posterior_moments.csv`, `M1_lactose_idata.nc` |
| 03 | `BruSafe_03_M2_soluble_partition.ipynb` | M2 | `VA_partition_gels.csv` if `USE_SYNTHETIC = False` | `M2_export.npz` |
| 04 | `BruSafe_04_M3_capsid_assembly.ipynb` | M3 | `M2_export.npz` (reference values are used if absent), `../data/derived_M3_dG_from_Tm.csv` | `M3_export.npz`, `M3_TEM_observable.npz` |
| 05 | `BruSafe_05_M5_absolute_quantification.ipynb` | M5 | inputs entered in the swap-point cells | none |

Section numbers in the notebooks refer to the CAPSID framework document, which is not part of this
repository.

## 2. Literature variants (`literature/`)

Each variant replaces the synthetic-data section of a production notebook with published
measurements. The data pack is read from `../../data/`.

| Notebook | Question addressed | Data |
|---|---|---|
| `BruSafe_01L_M1_literature_design.ipynb` | Does the ranking of IPTG, temperature and induction time predicted by M1 agree with a published factorial design? | `lit_M1_induction_design.csv` |
| `BruSafe_03L_M2_literature_arenin.ipynb` | What activation-energy parameter of the partition model is implied by the soluble fraction measured at two temperatures? | `lit_M2_partition_gels.csv` |
| `BruSafe_04L_M3_literature_dG.ipynb` | What contact energy follows from the published melting point of MS2 VLPs, and which M3 outputs depend on it? | `derived_M3_dG_from_Tm.csv`; `M2_export_literature.npz` from 03L |

The literature variants are older snapshots of the production notebooks in some respects. In
particular, notebook 01L keeps biomass in OD₆₀₀ units, whereas notebook 01 uses gDCW L⁻¹.

## 3. Exploratory notebooks (`exploratory/`)

Nothing here belongs to the production pipeline. The notebooks are earlier work on expression,
assembly and packaging, and they use a module numbering in which "M4" does not coincide with the
Module M4 of the framework. They are background for the M4 module, which the framework specifies
and which has no production notebook yet.

| Notebook | Language | Content |
|---|---|---|
| `Ejemplo_Capside_ensamblaje_N12.ipynb` | es | Kinetic master equation for a dodecahedral (N = 12) capsid after Zlotnick and Hagan; a short introductory example |
| `M3_expresion_colab.ipynb` | es | Expression of the coat protein to steady state; justification of D = P/2; temperature sweep (Colab-ready) |
| `M4_ensamblaje_colab.ipynb` | es | Explicit 90-step assembly cascade and the concentration sweep (Colab-ready) |
| `M4_L1_rate_equations.ipynb` | en | Layer 1: nucleation–elongation rate equations and the yield curve |
| `M4_L2_nucleation_phase_diagram.ipynb` | en | Layer 2: classical nucleation theory and the assembly phase diagram |
| `M4_L3_coarsegrained_MC.ipynb` | en | Layer 3: coarse-grained Monte Carlo of patchy subunits, which measures the contact ladder assumed in L1 |
| `M4_L4_packaging_signal_coassembly.ipynb` | en | Layer 4: stochastic co-assembly around cargo with packaging signals |
| `M4_bridge_CG_to_ODE_to_bench.ipynb` | en | Connects the CG bond census to the ODE and to benchtop observables (turbidity, TEM) |
| `M4_FIT_to_bench_data.ipynb` | en | Fits the contact energy, association rate and a recovery scale to bench data; profile likelihood |

Points to keep in mind when reading these notebooks:

- They use RT = 0.593 kcal mol⁻¹ (298 K), whereas the production notebooks use RT at 310.15 K
  (0.616 kcal mol⁻¹).
- The cascade of the exploratory notebooks predicts a yield that is peaked in the initial dimer
  concentration. The production cascade (notebook 04) uses a different contact schedule and the
  contact energy derived from the published melting point, and it produces a saturating trap and not a
  declining yield (§5 of notebook 04). The two sets of results should not be mixed.
- L1–L4, the bridge and the FIT notebook write their simulation source to disk in their first code
  cell. They refer to `M4_theory_compendium.pdf` and, for the FIT notebook, to
  `M4_bench_data_template.csv`. Neither file is included in this repository; without the template
  the FIT notebook falls back to a synthetic example.
- `Ejemplo_Capside_ensamblaje_N12.ipynb` is the only notebook with stored outputs in this folder.

## 4. Archived notebooks (`archive/`)

These are the direct predecessors of notebooks 01–05, under the module numbering used before the
framework document fixed M1–M5. Each starts with a banner cell that points to its successor; the
original text and outputs are otherwise unchanged. They are not to be run.

| Archived notebook | Successor |
|---|---|
| `BruSafe_M0_M4_bayesian_calibration_mM_longhorizon.ipynb` | `BruSafe_01_M1_expression_IPTG.ipynb` |
| `BruSafe_M0_M4_allolactose_24h_calibration.ipynb` | `BruSafe_02_M1_expression_lactose.ipynb` |
| `BruSafe_M2_5_soluble_partition.ipynb` | `BruSafe_03_M2_soluble_partition.ipynb` |
| `BruSafe_M3_capsid_assembly.ipynb` | `BruSafe_04_M3_capsid_assembly.ipynb` |
| `BruSafe_M6_keystone_CPP_and_scalability.ipynb` | `BruSafe_05_M5_absolute_quantification.ipynb` |

They still carry their original outputs, which is why the files are large. Several of them contain
results that the production notebooks have since corrected: the 300 min lactose horizon, the
uncensored Beta likelihood in the partition module, OD and band volume as primary observables, and a
contact energy of −4.2 kcal mol⁻¹ without a source.
