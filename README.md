# BruSafe · simulations and models

Computational work of the iGEM Tec-Chihuahua 2026 project on an MS2 virus-like-particle (VLP)
platform for mRNA delivery. The repository contains two arms:

- an **all-atom molecular dynamics (MD) campaign** with GROMACS, run on 10 workstations with one
  RTX 4070 each, and
- **CAPSID**, a chain of five kinetic and statistical models of the path from expression to purified
  VLP (`capsid/`).

Model deadline: 13 Oct 2026. Wiki freeze: 21 Oct 2026.

Wiki-ready text for both arms is in [`wiki/`](wiki/README.md).

## Repository layout

```
wiki/                  source pages for the wiki (overview, MD methods, CAPSID modules, data, limitations)

capsid/                the CAPSID modelling framework
  notebooks/             production, literature, exploratory and archived notebooks (catalogue in notebooks/README.md)
  data/                  literature data pack, bench-data templates and data dictionary
  docs/                  framework theory, results write-up, literature runs, M4 compendium
  figures/               figures produced by the notebooks, with an index
  results/logs/          cell-by-cell output of the reported runs
  tools/                 headless runner and calibration studies
  MANIFEST.csv           every notebook with group, module, size and title

config/systems.tsv     MD registry: what is simulated, how many replicas, and why
mdp/                   GROMACS parameter files (em, nvt, npt, npt_free, prod)
scripts/               prepare → equilibrate → produce → postprocess, plus run_night.sh and the GPU benchmark
analysis/              convergence checks, cross-system comparison, md_report.py (figures and wiki pages)
mmpbsa/                binding-energy (MM/PBSA) input and driver
docs/                  decision log, predictions, protocol (TUTORIAL), Windows setup, checklist, environment records
input/                 starting structures, one folder per system
```

Trajectories are excluded from version control. `.mdp`, `.top`, `.itp`, `.ndx` files, scripts and
analysis output are committed.

## Molecular dynamics

### Fixed parameters

| | |
|---|---|
| Force field | CHARMM36m (protein and nucleic acids) |
| Water | CHARMM-modified TIP3P |
| Salt | 150 mM KCl (`POT` / `CLA`) |
| Divalent ions | Mg²⁺ in RNA-containing systems only (see `docs/DECISIONS.md`) |
| Temperature | 310 K, V-rescale |
| Pressure | 1 bar, C-rescale, isotropic |
| Time step | 2 fs, constraints on bonds involving hydrogen |
| Box | rhombic dodecahedron, 1.2 nm padding |

### Platforms

Linux or WSL (preferred), macOS (CPU only, for preparation and analysis), and Windows Git Bash without
administrator rights (`docs/WINDOWS.md`).

### Quick start

```bash
# 0. record the state of this machine and commit the output
./scripts/00_env.sh

# 1. benchmark once, then fix the best mdrun flags for the whole campaign
./scripts/bench_gpu.sh runs/S1_wt_cc_apo/rep1/prod.tpr
export BRUSAFE_MDRUN_FLAGS="-nb gpu -pme gpu -bonded gpu -update gpu -ntmpi 1 -ntomp 8 -pin on"

# 2. plan the campaign across machines
python3 scripts/_registry.py plan --machines 10 --tier 1 --nsday <measured>

# 3. for each system and replica
./scripts/01_prepare.sh     S1_wt_cc_apo               # once per system
./scripts/02_equilibrate.sh S1_wt_cc_apo 1             # per replica
./scripts/03_production.sh  S1_wt_cc_apo 1 24          # restartable
./scripts/04_postprocess.sh S1_wt_cc_apo 1

# 4. the comparison that supports the main claim
./analysis/compare_systems.sh S1_wt_cc_apo S2_sccp_apo

# compact wiki pages: one per system, replicas pooled, one animation each
python3 analysis/md_report.py wiki --all
# full set of figures for one system, and cross-system overlays
python3 analysis/md_report.py report  S1_wt_cc_apo
python3 analysis/md_report.py compare S1_wt_cc_apo S2_sccp_apo

# 5. once the RNA-bound systems have been post-processed
PROT_GRP=<n> RNA_GRP=<n> ./mmpbsa/run_mmpbsa.sh S4_cp_pacdesign 1 4.0

# at any time: what is running
python3 scripts/_registry.py status --runs runs/
```

`03_production.sh` is restartable: after a crash or reboot, the identical command resumes from the
last checkpoint. The full protocol is in `docs/TUTORIAL.md`.

### Gates

1. **Equilibration gate** (`02_equilibrate.sh`): density ≈ 1000 kg m⁻³, temperature within 2 K of 310 K,
   and no remaining drift of the density.
2. **Duplicate-run guard** (`guard_existing`): a finished replica is not overwritten.
3. **Decomposition gate** (`run_mmpbsa.sh`): the per-residue hot spots must recover the known
   operator-contact residues before any binding energy is interpreted.

### Reporting rules

- **ΔΔG only.** Absolute MM/PBSA binding free energies are commonly off by 5–20 kcal mol⁻¹, and more
  so for charged protein–RNA interfaces.
- **ΔΔG = −RT ln(ratio)**, with RT = 0.616 kcal mol⁻¹ at 310 K. One order of magnitude in K<sub>d</sub>
  corresponds to about 1.4 kcal mol⁻¹.
- **Error bars come from the spread among replicas**, not from the number of frames. The naive
  standard error over correlated frames underestimates the uncertainty by roughly an order of magnitude
  (`convergence.py block` prints the factor for the data at hand).
- **State whether entropy is included.** Without −TΔS the quantity is an effective interaction
  enthalpy and not a ΔG.
- **Register predictions before the readouts** in `docs/PREDICTIONS.md`.

## CAPSID

Five modules (expression, soluble partition, *T* = 3 assembly, cargo packaging, absolute
quantification) are fitted each on its own observable, and only a concentration is passed from one to
the next. Module M4 (packaging) is specified but not yet implemented as a production notebook. See
[`capsid/README.md`](capsid/README.md) for the notebooks, how to run them, and the known limitations.
