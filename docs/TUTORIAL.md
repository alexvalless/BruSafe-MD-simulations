# Running the campaign — protocol

The full workflow should be run once from start to finish on **S1** before any other system is
started. If S1 works, every other system needs the same four commands with a different name.

Terminology: a *system* is a molecular setup (`S1_wt_cc_apo`). A *replica* is one independent run
of that system (`rep1`, `rep2`, `rep3`).

---

## 0. What the scripts do not do

The scripts handle solvation, ionisation, equilibration, production, PBC correction, analysis and
bookkeeping.

They do **not** build the starting structures. Selecting the C/C chain pair from 2MS2, modelling
the single-chain linker and mutating the RNA to the designed pac sequence are scientific decisions,
they are where the errors that matter arise, and they are deliberately left to the user. Section 2
describes them.

---

## 1. One-time setup on every machine

### 1.1 GROMACS

```bash
gmx --version | head -3
```

Version **2023 or newer** is required: `gmx dssp` does not exist in earlier releases, and
GPU-resident update performs much better from 2023 on. The build must be CUDA-enabled:

```bash
gmx --version | grep -i 'GPU support'      # expected: CUDA
```

**Every machine must report the identical version.** Replicas of one system run with different
GROMACS versions are not comparable, which silently invalidates the error bars. If one machine has
2023.3 and the others 2024.2, this must be corrected before anything is run.

### 1.2 CHARMM36m force field

GROMACS does not include it. The GROMACS port is downloaded from the MacKerell laboratory, unpacked,
and `GMXLIB` is set to its parent directory:

```bash
mkdir -p ~/ff && cd ~/ff
# unpack so that ~/ff/charmm36-jul2022.ff/ exists
echo 'export GMXLIB=$HOME/ff' >> ~/.bashrc && source ~/.bashrc
ls $GMXLIB/charmm36*.ff/forcefield.doc     # must exist
```

If the directory name differs from `charmm36-jul2022`, the `-ff` argument in `scripts/01_prepare.sh`
has to be updated, and the change recorded in `docs/DECISIONS.md`.

### 1.3 gmx_MMPBSA (one or two machines only)

It is not needed for production and is installed only where M3.5 will be run:

```bash
conda create -n mmpbsa -c conda-forge python=3.10 ambertools mpi4py -y
conda activate mmpbsa
pip install gmx_MMPBSA
gmx_MMPBSA --version
```

### 1.4 Repository and environment record

```bash
git clone <remote> brusafe-md && cd brusafe-md
chmod +x scripts/*.sh analysis/*.sh mmpbsa/*.sh
./scripts/00_env.sh
```

This writes `docs/env_<hostname>.txt`. The file should be committed from every machine: it
documents the computing environment of the campaign, and it is where a disagreement between two
machines shows up first.

### 1.5 Location of the trajectories

By default trajectories are written to `runs/` inside the repository. With shared storage, every
machine can be pointed at it so that the duplicate-run guard works across machines:

```bash
echo 'export BRUSAFE_RUNS=/mnt/shared/brusafe/runs' >> ~/.bashrc
```

Without shared storage the guard protects only against duplicates on a single machine, and the
assignment of replicas to machines has to be tracked separately.

---

## 2. Input structures (the manual part)

The step-by-step version, including the CHARMM-GUI settings and the `scripts/00_build_inputs.py` helper that assembles the dimers and mutates the RNA, is in `docs/TOPOLOGIES.md`. This section gives the reasoning behind it.

Each system expects one of the following:

```
input/<system>/<system>.pdb        # pdb2gmx route  (apo systems)
input/<system>/gromacs/            # CHARMM-GUI Solution Builder output (RNA systems)
```

### S1 — WT C/C dimer from 2MS2

The asymmetric unit contains three quasi-equivalent chains and therefore **two physically different
dimers**: the symmetric C/C pair and the asymmetric A/B pair. The free coat-protein dimer in
solution resembles the C/C conformer, whereas A/B is the conformer induced by assembly and RNA
binding. Simulating the wrong one gives a baseline that is already locked into the assembly state.

Chain identifiers should not be taken at face value. The biological assembly is inspected to
identify the pair in which **both** FG loops are in the extended conformation. The dimer is then
extracted:

```bash
mkdir -p input/S1_wt_cc_apo
pdb_selchain -C,C 1msc.pdb | pdb_delhetatm | pdb_tidy > input/S1_wt_cc_apo/S1_wt_cc_apo.pdb
# or without pdb-tools:
#   grep '^ATOM' 1msc.pdb | awk 'substr($0,22,1)=="C"' > ...
```

Conditions for accepting a PDB file as input:

- [ ] crystallographic waters removed
- [ ] alternate locations resolved to a single conformer
- [ ] cryoprotectants and buffer components removed
- [ ] no missing residues in the FG loop (built if absent)
- [ ] termini chosen deliberately (`01_prepare.sh` runs `pdb2gmx -ter` interactively)

### S2 — single-chain dimer

The model is built with AlphaFold3 or by grafting the linker onto the 2MS2 dimer. **Before any
simulation**, the model is superposed onto the 2MS2 C/C dimer and the backbone RMSD is recorded:

```bash
gmx confrms -f1 input/S1_wt_cc_apo/S1_wt_cc_apo.pdb \
            -f2 input/S2_sccp_apo/S2_sccp_apo.pdb -o fit.pdb
```

This value belongs on the wiki in any case. If the two CP subunits are in a non-native relative
orientation, every later S1-versus-S2 comparison measures an AlphaFold artefact and not an effect of
the linker, and no amount of sampling corrects that.

### S3, S4, S8, S9 — RNA-bound systems

The CP–operator **cocrystal** is used. An RNA hairpin should not be docked de novo when a solved
complex exists. The designed pac and the calibration variants are built by mutating the cocrystal
RNA in place.

For anything containing RNA, CHARMM-GUI Solution Builder is preferred over the `pdb2gmx` route:

```bash
./scripts/01_prepare.sh S4_cp_pacdesign --source charmm-gui
```

There are two reasons. The `pdb2gmx` residue mapping for RNA under the charmm36 port is fragile,
and CHARMM-GUI handles Mg²⁺ placement properly, which matters more than it may seem (see §7).

---

## 3. S1 from start to finish

```bash
./scripts/01_prepare.sh     S1_wt_cc_apo --replica 1   # once per replica (own ion placement)
./scripts/02_equilibrate.sh S1_wt_cc_apo 1    # once per replica
```

`01_prepare.sh` prints the index groups it has built. `Backbone`, `C-alpha`, `SOLU` and `SOLV`
should be present; the script stops if any is missing, because `04_postprocess.sh` requests them by
name and an early failure is much less costly than one several days later.

`02_equilibrate.sh` runs EM → NVT 200 ps restrained → NPT 1 ns restrained → NPT 2 ns free, and then
applies the **equilibration gate**: density close to 1000 kg/m³, temperature within 2 K of 310, and
density no longer drifting. It does not exit successfully otherwise. A drifting box indicates an
unsettled system whatever the mean value is, and the gate detects this where a mean alone would not.

About 20 minutes should be spent here inspecting `density.xvg` and the RMSD, since this is the last
inexpensive point at which a problem can be noticed.

---

## 4. Benchmark once, then lock the configuration

A `prod.tpr` is needed for the benchmark, so one is generated by hand from the equilibrated replica:

```bash
cd runs/S1_wt_cc_apo/rep1
gmx grompp -f ../../../mdp/prod.mdp -c npt_free.gro -t npt_free.cpt \
           -p topol.top -n index.ndx -o prod.tpr
cd -

./scripts/bench_gpu.sh runs/S1_wt_cc_apo/rep1/prod.tpr
rm runs/S1_wt_cc_apo/rep1/prod.tpr    # 03_production.sh regenerates it from the registry
```

The `rm` is necessary. `03_production.sh` sets `nsteps` from the `ns` column of
`config/systems.tsv`, and the RNA systems run 30 ns instead of 100; a hand-made tpr would silently override
that.

The sweep tests three offload configurations at several thread counts and then runs two concurrent
jobs on one GPU. At ~35k atoms a single run usually does not saturate an RTX 4070, so the
**aggregate throughput of two concurrent runs is often higher than that of one run**, and this is
the number that determines the schedule.

The chosen configuration is then fixed on every machine:

```bash
echo 'export BRUSAFE_MDRUN_FLAGS="-nb gpu -pme gpu -bonded gpu -update gpu -ntmpi 1 -ntomp 8 -pin on"' >> ~/.bashrc
```

The selected flags and the date are recorded in `docs/DECISIONS.md`. Every replica in the campaign
must use identical flags, otherwise the timings are not comparable.

If one machine benchmarks 25% below the others, thermal or power throttling should be excluded
before the software is suspected:

```bash
nvidia-smi -q -d PERFORMANCE | grep -A5 'Clocks Event Reasons'
```

Poor airflow in OEM desktops is a real and common source of spread among nominally identical
machines.

---

## 5. Planning across the ten machines

```bash
python3 scripts/_registry.py plan --machines 10 --tier 1 --nsday <measured value>
```

The output is a cost table and a list of jobs per machine. Tier 1 amounts to about 1,400 ns-equivalent including equilibration; at 400 ns/day per GPU
(hydrogen mass repartitioning, 4 fs) this is well under one day on ten GPUs. Compute is therefore not
the limiting factor, which is the reason for running this early.

The job lists are saved and one is assigned to each machine:

```bash
python3 scripts/_registry.py plan --machines 10 --tier 1 --nsday <measured value> > docs/PLAN.txt
git add docs/PLAN.txt && git commit -m "campaign plan, tier 1"
```

---

## 6. The production loop

For each replica, on the assigned machine:

```bash
./scripts/01_prepare.sh     S1_wt_cc_apo --replica 2   # skipped if build_rep2/ already exists
./scripts/02_equilibrate.sh S1_wt_cc_apo 2
./scripts/03_production.sh  S1_wt_cc_apo 2 23    # 23 = wall-clock hours per chunk
./scripts/04_postprocess.sh S1_wt_cc_apo 2
```

`03_production.sh` is **restartable**. When it reaches the wall-clock limit or the machine reboots,
the identical command resumes from the last checkpoint. A simple loop runs a replica to completion
overnight:

```bash
until ./scripts/03_production.sh S1_wt_cc_apo 2 23; do sleep 10; done
```

The script also writes `PROVENANCE.txt` next to the trajectory, containing the host, GROMACS
version, mdrun flags, seed, and start and finish times. This file makes a run traceable weeks later,
when nobody remembers which machine did what.

Overall progress can be checked at any time:

```bash
python3 scripts/_registry.py status --runs runs/
```

The check should be made **daily**. With ten unmanaged desktops the realistic failure is not a lack
of compute but two machines running the same replica while another system is never launched.

---

## 7. Failure modes

**Equilibration gate fails.** `density.xvg` is plotted. A slow, monotonic drift usually means that
the box has not settled, and `npt_free` is extended and rerun. A sudden jump or a crash when the
restraints are released usually indicates a poor starting structure, and the work returns to §2. The
gate should not be bypassed: production from an unequilibrated box wastes GPU-days and the problem
only becomes visible at the analysis stage.

**`already holds a completed run`.** The duplicate guard has worked as intended. After confirming
that nobody else is running that replica, `BRUSAFE_FORCE=1 ./scripts/03_...` overrides it if
overwriting is really intended.

**`gmx dssp` unavailable.** GROMACS is older than 2023. Either upgrade, or use `do_dssp` with an
external DSSP binary. Everything else still runs.

**`no 'RNA' group`.** `make_ndx` did not recognise the nucleic-acid residue names. The group is
added manually:

```bash
gmx make_ndx -f runs/<sys>/build/solv_ions.gro \
             -n runs/<sys>/build/index.ndx -o runs/<sys>/build/index.ndx
```

**Mg²⁺ in the interface.** `04_postprocess.sh` writes `mg_interface.txt` for RNA systems (the fraction
of frames in which each ion is within 0.5 nm of the RNA, of the protein, and of both; a non-zero value
for both indicates that the ion sits in the interface) together with the raw distance files
`mindist_mg_rna.xvg` and `mindist_mg_protein.xvg`. These should be inspected. Water exchange at Mg²⁺
takes microseconds, so on a 30–100 ns trajectory the ions
remain essentially where they were placed; random `genion` placement adds noise and no realism. An
ion located in the protein–RNA interface contaminates every M3.5 energy from that replica. The
preferred placement order is crystallographic sites, then CHARMM-GUI, then `genion` as a documented
last resort. This is also a suitable question for the external MD advisor and worth documenting for
the wiki.

---

## 8. Weekly rhythm

- **Daily:** `_registry.py status`. Stalled runs are restarted the same day.
- **Twice a week:** analysis output, `docs/env_*.txt` and `PROVENANCE.txt` are committed. `.xtc`
  files are never committed (they are gitignored); trajectories are backed up separately, since
  their loss is the one unrecoverable failure.
- **27 September:** cut-off for new systems. Systems not launched by then are out of scope.
  Compute was never the limiting factor.

---

## 9. The comparison that supports the main claim

The RMSD of a single system demonstrates little. The claim rests on the comparison:

```bash
./analysis/compare_systems.sh S1_wt_cc_apo S2_sccp_apo
```

This computes the **essential-subspace overlap** between the principal motions of the two systems.
Over the first 10 eigenvectors, a value above 0.7 supports the statement that the linker does not
alter the dynamics of the dimer, and a value below 0.5 indicates a real change. Two RMSF curves that
look similar by eye are not evidence; the overlap is.

The script also prints the inter-replica RMSF agreement within each system, and **this should be
read first**. A difference between S1 and S2 is meaningful only if it exceeds the spread among
replicas of the same system. Many over-interpreted MD results fail at exactly this point.

One practical point: S2 contains additional linker residues, so the two systems do not have the same
number of atoms. The covariance analysis should be restricted to the CP backbone atoms present in
both before the comparison is made.

---

## 10. M3.5

M3.5 is run only after `04_postprocess.sh` has been run on the RNA-bound system. The group numbers
are identified first:

```bash
gmx make_ndx -f runs/S4_cp_pacdesign/build/solv_ions.gro -o /dev/null  # note the numbers, then q

export PROT_GRP=<n> RNA_GRP=<n>
./mmpbsa/run_mmpbsa.sh S4_cp_pacdesign 1 4.0
```

The output is checked **in this order**, which is fixed:

1. `FINAL_DECOMP.dat`: do the per-residue hot spots recover the known operator-contact residues of
   MS2 CP? If not, the setup is wrong and no energy computed from it is meaningful. This check takes
   a day and gates everything that follows.
2. The calculation is repeated at `indi` = 1, 2 and 4, and the spread is reported as part of the
   error bar. A conclusion that changes between 2 and 4 is not robust.
3. The calculation is repeated on two disjoint frame windows.
4. Only the **ΔΔG between systems** is reported. Absolute MM/PBSA binding free energies are
   routinely off by 5–20 kcal/mol, and more so for charged protein–RNA interfaces; they are not
   quoted.
5. ΔΔG = −RT ln(ratio), RT = 0.616 kcal/mol at 310 K. One order of magnitude in K_d corresponds to
   about 1.4 kcal/mol.

---

## 10b. Figures and animations for the wiki

`analysis/md_report.py` turns the output of `04_postprocess.sh` into figures. It requires a few
additional Python packages, on the analysis machine only:

```bash
pip install numpy matplotlib plotly MDAnalysis   # MDAnalysis only for animations
```

```bash
python3 analysis/md_report.py report  S1_wt_cc_apo              # all replicas overlaid
python3 analysis/md_report.py compare S1_wt_cc_apo S2_sccp_apo  # cross-system overlays
python3 analysis/md_report.py animate --top runs/S1_wt_cc_apo/rep1/analysis/clean.gro \
        --traj runs/S1_wt_cc_apo/rep1/analysis/clean.xtc -o s1_rep1.html
```

**For the wiki, the `wiki` command is the starting point.** It writes one compact page per system to
`analysis_out/wiki/<system>/html/index.html` (summary table; RMSD, RMSF, Rg and SASA with all
replicas overlaid and mean ± SD; RNA RMSD and Mg²⁺ distance for RNA systems; one animation, with
`--anim-rep N` selecting the replica) and an `analysis_out/wiki/index.html` that links them:

```bash
python3 analysis/md_report.py wiki --all                    # every system with analysis output
python3 analysis/md_report.py wiki S4_cp_pacdesign S5_cp_scramble
```

Systems with nothing to plot are skipped with a message. `report` produces the complete set, for
the supplement or for checking convergence.

`report` writes to `analysis_out/<system>/`:

| Output | Content |
|---|---|
| `figures/*.png`, `*.svg` | RMSD, RNA RMSD, RMSF (replicas + mean ± SD), Rg, SASA, H-bonds, Mg²⁺ distance, DSSP map and counts, PCA (variance, PC1 with cosine content, PC1×PC2), block-averaging SEM, distributions |
| `html/index.html` | dashboard: summary table and every plot, interactive |
| `html/<plot>.html` | one plot per page, for an `<iframe>` on the wiki |
| `html/traj_rep<N>.html` | 3Dmol.js trajectory player: play/pause, slider, cartoon/sticks, colouring by chain, N→C, secondary structure or RMSF |
| `summary.csv`, `summary.md` | per-replica means and mean ± SD across replicas |

**For the iGEM wiki**, the whole `html/` folder is uploaded. The default `--js local` places
`plotly.min.js` and `3Dmol-min.js` next to the pages, because the wiki does not allow external CDNs
(3Dmol.js is vendored in `analysis/vendor/`, BSD-3). `--js inline` instead produces self-contained
single files.

Useful options: `--frames 100` (frames per movie; fewer frames give a smaller file), `--smooth 5`
(averages neighbouring frames to damp thermal jitter in the movie; this should be stated in the
caption), `--sel "name CA P"` (Cα/P atoms only, for large systems such as S10), `--equil-ns 50`
(relaxation window excluded from all statistics), and `--no-anim`.

The animation must be made from `clean.xtc`: a trajectory that is not PBC-whole gives a movie with
bonds stretched across the box, and the script warns about this.

---

## 11. Before the wiki freeze

- [ ] `docs/env_*.txt` committed from every machine
- [ ] `docs/DECISIONS.md` complete, including the locked mdrun flags and any deviation from the
      defaults
- [ ] `docs/PREDICTIONS.md` timestamped **before** the dPCR readout
- [ ] Systems table: atom counts, box, lengths, replicas, achieved ns/day
- [ ] Block-averaging factor from the campaign's own data quoted in the methods section
      (`convergence.py block` prints it); this anticipates the most common criticism of
      student MD analyses
- [ ] An explicit section on what was not validated; stating the limitations in the authors' own
      words is preferable to leaving them to a reviewer
