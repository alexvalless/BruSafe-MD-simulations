# Reproducibility

## CAPSID notebooks

**Environment.** Python 3.11 with NumPy ≥ 2.0, SciPy ≥ 1.14, pandas, matplotlib, numba,
PyMC ≥ 5.20, ArviZ and h5netcdf (`capsid/environment.yml`; `capsid/requirements.txt` for pip).

```bash
conda env create -f capsid/environment.yml && conda activate capsid
cd capsid/notebooks && jupyter lab
```

**Working directory.** Each notebook is run from its own folder. The notebooks exchange files
through the working directory (for example `M2_export.npz`, from notebook 03 to notebook 04) and the
literature variants read the data pack through relative paths. Notebooks within one folder are run in
numerical order.

**Cost.** `FAST_MODE = True` gives a usable posterior in a few minutes. The numbers in the pages
come from `FAST_MODE = False`, with 2000 draws × 4 chains (or 2000 SMC particles).

| Notebook | Runtime (full settings) |
|---|---|
| 01 · M1, IPTG | ~16 min |
| 02 · M1, lactose | ~14 min |
| 03 · M2 | ~25 min |
| 04 · M3 | ~15 s |
| 05 · M5 | ~1 s |

**Random seeds.** The seeds are fixed in the notebooks (for example 2026 in notebook 01), so a rerun
with the same environment reproduces the numbers.

**Logs.** The cell-by-cell output of the runs that produced the reported numbers is in
`capsid/results/logs/`.

**Calibration studies.** `capsid/tools/sbc_1d.py` (exact one-dimensional posteriors on a grid),
`sbc_m2.py` (joint four-dimensional posterior, with sampler diagnostics) and `diag_m2.py` (a long
reference chain against the SMC posterior). The numbers quoted in the text came from 64 calibration
replicates; the rank histogram shown in the figure comes from the 16 replicates of the notebook default.

**Headless execution.**

```bash
cd capsid/notebooks
python ../tools/exec_nb.py BruSafe_01_M1_expression_IPTG.ipynb run01
```

## Molecular dynamics

**Software.** GROMACS 2023 or later with CUDA support (identical on all machines), the CHARMM36m port
(`charmm36-jul2022.ff`), and gmx_MMPBSA with AmberTools on one or two analysis machines. Linux and
WSL are preferred; macOS (analysis only) and Windows Git Bash are supported (`docs/WINDOWS.md`).

**Record of the environment.** `scripts/00_env.sh` writes `docs/env_<hostname>.txt` with the hardware,
software versions, force-field location and locked mdrun flags of each machine, and the files are
committed. Each production replica also writes a `PROVENANCE.txt` (host, GROMACS version, mdrun
flags, seed, start and end times).

**Determinism.** Velocity seeds are a deterministic hash of `<system>_rep<N>`, so every equilibration
can be reproduced from the repository.

**Workflow.**

```bash
./scripts/00_env.sh
./scripts/bench_gpu.sh runs/S1_wt_cc_apo/rep1/prod.tpr      # once; the result fixes the mdrun flags
./scripts/01_prepare.sh     S1_wt_cc_apo
./scripts/02_equilibrate.sh S1_wt_cc_apo 1
./scripts/03_production.sh  S1_wt_cc_apo 1 24
./scripts/04_postprocess.sh S1_wt_cc_apo 1
./analysis/compare_systems.sh S1_wt_cc_apo S2_sccp_apo
python3 analysis/md_report.py wiki --all
```

Trajectories, checkpoints and run-time files are excluded from version control (`.gitignore`) and
are backed up separately. Topologies, index files, `.mdp` files, scripts and analysis output are
committed.
