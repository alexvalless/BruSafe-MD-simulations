# BruSafe · CAPSID

**C**omputational **A**ssembly **P**latform for **S**tructural **I**mmunogenicity and **D**elivery.
Modelling framework of the iGEM Tec-Chihuahua 2026 project: an MS2 virus-like particle (VLP)
that delivers mRNA.

The framework consists of five modules chained so that the output of one module is the input of the
next: expression → soluble partition → *T* = 3 assembly → cargo packaging → absolute
quantification. Each module is fitted to its own observable, and the only quantity passed to the
next module is a concentration.

![Model overview](figures/fig28_wiki_model_overview.png)

## Contents

```
notebooks/              the five production modules, in pipeline order
notebooks/literature/   the same modules confronted with published measurements
notebooks/exploratory/  earlier assembly and packaging work (separate numbering)
notebooks/archive/      predecessors of the production notebooks, kept for provenance
data/                   literature data pack (CSV) used by the notebooks, with a data dictionary
figures/                the figures produced by the notebooks (150 dpi), with an index
results/logs/           cell-by-cell output of the runs reported in the write-up
tools/                  headless notebook runner and the calibration studies
MANIFEST.csv            every notebook, with group, module, language, size and title
environment.yml         conda environment (requirements.txt for pip)
```

A catalogue of the notebooks with their inputs and outputs is in
[`notebooks/README.md`](notebooks/README.md).

## Production notebooks

The notebooks are run in order. Notebooks 03 and 04 are linked through the file `M2_export.npz`. The other
links between modules are made through values that are entered in the receiving notebook: notebook 03 uses
a reproduction of the M1 forward model at fixed parameters and does not read the posterior draws of
notebooks 01–02, and the process inputs of notebook 05 are entered in its swap-point cell.

| Notebook | Module | Content | Runtime |
|---|---|---|---|
| `BruSafe_01_M1_expression_IPTG.ipynb` | M1 | inducer uptake, *lac* de-repression and CP synthesis; SMC fit to a biomass and protein time course under IPTG | ~16 min |
| `BruSafe_02_M1_expression_lactose.ipynb` | M1 | the same system under lactose autoinduction, with glucose, LacY and β-galactosidase explicit | ~14 min |
| `BruSafe_03_M2_soluble_partition.ipynb` | M2 | partition between folding and inclusion bodies as a function of synthesis flux and temperature, including the calibration study | ~25 min |
| `BruSafe_04_M3_capsid_assembly.ipynb` | M3 | 90-step nucleation–elongation cascade, critical concentration, feed versus bolus, TEM observable | ~15 s |
| `BruSafe_05_M5_absolute_quantification.ipynb` | M5 | propagation from culture to vial: harvest titre, cargo copies per particle, error budget | ~1 s |

Module M4 (CP–*pac* energetics and packaging calibration) is specified in the framework document
but has no production notebook yet. The notebooks in `notebooks/exploratory/` are earlier work
towards it and use a different module numbering; they are background and not part of the pipeline.

## Literature variants

Each variant is the corresponding production notebook with its synthetic-data section replaced by
published measurements. The remaining sections are unchanged, so differences in the output are due
to the data.

| Notebook | Data | Source |
|---|---|---|
| `BruSafe_01L_M1_literature_design.ipynb` | 15-run Box–Behnken design, three proteins in BL21(DE3) | doi:10.1038/s41598-019-53200-7 |
| `BruSafe_03L_M2_literature_arenin.ipynb` | soluble and insoluble protein at two induction temperatures | doi:10.1186/s40643-025-00986-2 |
| `BruSafe_04L_M3_literature_dG.ipynb` | thermal denaturation midpoints of MS2 VLPs | doi:10.1021/acs.biochem.4c00225 |

## Running the notebooks

```bash
conda env create -f environment.yml && conda activate capsid
cd notebooks && jupyter lab        # interactive
```

The notebooks exchange files (`M1_*_posterior_draws.csv`, `M2_export.npz`, …) through the
**working directory**, and the literature variants read the data pack through relative paths
(`../data/` from `notebooks/`, `../../data/` from `notebooks/literature/`). Each notebook should
therefore be run with its own folder as the working directory, which is the default in Jupyter, and
the notebooks of one folder should be run in numerical order.

For unattended runs:

```bash
cd notebooks
python ../tools/exec_nb.py BruSafe_01_M1_expression_IPTG.ipynb run01
```

`tools/exec_nb.py` executes a notebook in-process instead of through a Jupyter kernel, saves every
figure to `figs/`, and logs the timing and output of each cell. It was written because the
environment in which the notebooks were developed could not start a kernel; on a normal machine
Jupyter is sufficient, and the runner is useful only for unattended runs.

Each production notebook has one `USE_SYNTHETIC = True` flag and one marked data cell. Setting the
flag to `False` makes the notebook read a CSV file with the columns listed in that cell; nothing
else needs to change when real data become available.

Two settings control the cost. `FAST_MODE = True` gives a usable posterior in a few minutes; the
numbers in the write-up come from `FAST_MODE = False` (2000 draws × 4 chains, or 2000 SMC
particles). `R_SBC` in notebook 03 sets the number of calibration replicates: 16 runs take about two
and a half minutes, and the reported figure used 64, through `tools/sbc_1d.py`.

## Observables

The modules use physical quantities and not instrument readings:

- **biomass** in gDCW L⁻¹, from OD₆₀₀ × 0.40 or from wet pellet weight × 0.25
- **coat protein** in mg L⁻¹, by densitometry against a purified standard curve with one slope per gel
- **soluble and insoluble fractions** as two concentrations in mg L⁻¹, not as a band ratio
- **correctly formed fraction** from TEM: of at least 200 particles, those with a diameter of
  27 ± 2 nm, with a binomial interval
- **particle concentration** from SEC with a protein assay, and **cargo copies** by dPCR

## Tools

| File | Purpose |
|---|---|
| `tools/exec_nb.py` | headless notebook execution with figure capture and per-cell logging |
| `tools/sbc_1d.py` | simulation-based calibration for M2 with exact one-dimensional posteriors on a grid, without MCMC |
| `tools/sbc_m2.py` | the same check over the joint four-dimensional posterior, with sampler diagnostics for each replicate |
| `tools/diag_m2.py` | long reference Metropolis–Hastings chain compared with the SMC posterior, to separate sampler error from model error |

The calibration scripts read `BruSafe_03_M2_soluble_partition.ipynb` from the working directory and write
to `handoff/`, so they are run from `notebooks/` after `mkdir -p handoff`
(for example `python ../tools/sbc_1d.py`).

The calibration study was not an optional extra. It identified a real defect in the observation
model of M2: gels reading a soluble fraction of exactly 1 were treated as infinitely informative
instead of as censored observations, and a conventional recovery check had accepted this for
weeks.

## Known limitations

- The contact schedule in M3 is a placeholder with the correct total. It moves the kinetic capsid
  fraction at 18 h between 0.003 and 0.930, and it is the largest unknown in the framework.
  Enumerating the contacts from PDB 2MS2 and 1AQ3 requires no bench work.
- Held-out dose coverage is 50% (IPTG) and 60% (lactose) against a nominal 94%: the model
  interpolates among the fitted doses and is overconfident at an unseen dose. A Kennedy–O'Hagan
  discrepancy term is required before any forecast interval at an unseen condition is quoted.
- The recovery efficiencies of M5 are assumed and not measured. They contribute 55% of the variance
  in cargo copies per particle, and the three controls that would measure them are the most valuable
  experiment in the plan.
- The theory compendium for the exploratory M4 notebooks (`M4_theory_compendium.pdf`), the bench-data
  template (`M4_bench_data_template.csv`) and the framework document to which the section numbers
  refer are not part of this repository.

## Licence and citation

The code is released for the iGEM 2026 competition by Tec-Chihuahua. If it is adapted, the module
structure and the identifiability audits transfer to any VLP platform; the geometry constants
(*N* = 90, *C*<sub>tot</sub> = 180) and the priors do not.
