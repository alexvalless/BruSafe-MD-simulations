# BruSafe MD

Molecular dynamics campaign for the iGEM 2026 BruSafe MS2 VLP platform.
GROMACS on 10 × RTX 4070. Model deadline **13 Oct 2026**, wiki freeze **21 Oct**.

## Locked parameters

| | |
|---|---|
| Force field | CHARMM36m (protein + nucleic acid) |
| Water | CHARMM-modified TIP3P |
| Salt | 150 mM KCl (`POT` / `CLA`) |
| Divalent | Mg²⁺, RNA-containing systems only — see the caveat below |
| Temperature | 310 K, V-rescale |
| Pressure | 1 bar, C-rescale, isotropic |
| Timestep | 2 fs, h-bond constraints |
| Box | rhombic dodecahedron, 1.2 nm padding |


## Layout

```
config/systems.tsv     the registry: what runs, how many replicas, and why
mdp/                   em, nvt, npt, npt_free, prod 
scripts/               prepare → equilibrate → produce → postprocess
analysis/              convergence checks and cross-system comparison
mmpbsa/                Binding energetics
docs/                  decision log, pre-registered predictions, env dumps
```

Trajectories are gitignored. Commit `.mdp`, `.top`, `.itp`, `.ndx`, scripts and
analysis output. 

## Quickstart

```bash
# 0. record this machine's exact state; commit the output
./scripts/00_env.sh

# 1. benchmark ONCE, then lock the winning mdrun flags for the whole campaign
./scripts/bench_gpu.sh runs/S1_wt_cc_apo/rep1/prod.tpr
export BRUSAFE_MDRUN_FLAGS="-nb gpu -pme gpu -bonded gpu -update gpu -ntmpi 1 -ntomp 8 -pin on"

# 2. plan the campaign across machines
python3 scripts/_registry.py plan --machines 10 --tier 1 --nsday <measured>

# 3. per system, per replica
./scripts/01_prepare.sh     S1_wt_cc_apo               # once per system
./scripts/02_equilibrate.sh S1_wt_cc_apo 1             # per replica
./scripts/03_production.sh  S1_wt_cc_apo 1 24          # restartable
./scripts/04_postprocess.sh S1_wt_cc_apo 1

# 4. the comparison that carries the actual claim
./analysis/compare_systems.sh S1_wt_cc_apo S2_sccp_apo

# 5. Once RNA-bound systems are postprocessed
PROT_GRP=<n> RNA_GRP=<n> ./mmpbsa/run_mmpbsa.sh S4_cp_pacdesign 1 4.0

# anytime: what is actually running
python3 scripts/_registry.py status --runs runs/
```

`03_production.sh` rerun the identical command after a crash
or reboot and resumes from the last checkpoint.

## Gates

1. **Equilibration gate** (`02_equilibrate.sh`) — density ≈ 1000 kg/m³,
   temperature within 2 K of 310, and density no longer drifting.
2. **Duplicate-run guard** (`guard_existing`) — refuses to overwrite a
   finished replica. 
3. **Decomposition gate** (`run_mmpbsa.sh`) — per-residue hot spots must
   recover the known operator-contact residues before any binding energy is
   interpreted.

## Reporting rules

- **ΔΔG only.** Absolute MM/PBSA binding free energies are routinely off by
  5–20 kcal/mol and worse for charged protein–RNA interfaces. 
- **ΔΔG = −RT ln(ratio)**, RT = 0.616 kcal/mol at 310 K. One order of
  magnitude in K_d ≈ 1.4 kcal/mol.
- **Error bars come from replica spread**, not from frame count. The naive SEM
  over correlated frames underestimates by roughly an order of magnitude
  (`convergence.py block` prints the factor for your own data) and is the
  clearest tell of an inexperienced analysis.
- **State whether entropy is included.** Without −TΔS it is an effective
  interaction enthalpy, not a ΔG.
- **Register predictions before readouts** in `docs/PREDICTIONS.md`.
