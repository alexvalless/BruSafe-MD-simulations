# Running the campaign — tutorial

Follow this once end-to-end on **S1** before touching any other system. If S1
works, everything else is the same four commands with a different name.

Terminology: a *system* is a molecular setup (`S1_wt_cc_apo`). A *replica* is
one independent run of that system (`rep1`, `rep2`, `rep3`).

---

## 0. What the scripts do not do

Be clear about the boundary before you start. The scripts handle solvation,
ionisation, equilibration, production, PBC cleanup, analysis and bookkeeping.

They do **not** build your starting structures. Choosing the C/C chain pair out
of 2MS2, modelling the single-chain linker, mutating the RNA to your designed
pac sequence — that is science, it is where the mistakes that matter live, and
it is deliberately left to you. Section 2 covers it.

---

## 1. One-time setup, on every machine

### 1.1 GROMACS

```bash
gmx --version | head -3
```

Needs **2023 or newer** — `gmx dssp` does not exist before that, and GPU-resident
update is much better from 2023 on. Build must be CUDA-enabled:

```bash
gmx --version | grep -i 'GPU support'      # expect: CUDA
```

**Every machine must report the identical version.** Mixed GROMACS versions
across replicas of the same system makes them non-comparable and quietly
invalidates your error bars. If one box has 2023.3 and the rest have 2024.2,
fix it before running anything.

### 1.2 CHARMM36m force field

GROMACS does not ship it. Download the GROMACS port from the MacKerell lab,
unpack, and point `GMXLIB` at the parent directory:

```bash
mkdir -p ~/ff && cd ~/ff
# unpack so that ~/ff/charmm36-jul2022.ff/ exists
echo 'export GMXLIB=$HOME/ff' >> ~/.bashrc && source ~/.bashrc
ls $GMXLIB/charmm36*.ff/forcefield.doc     # must exist
```

If the directory name differs from `charmm36-jul2022`, update the `-ff` argument
in `scripts/01_prepare.sh` **and** record the change in `docs/DECISIONS.md`.

### 1.3 gmx_MMPBSA — one or two machines only

Not needed for production. Install it where you will run M3.5:

```bash
conda create -n mmpbsa -c conda-forge python=3.10 ambertools mpi4py -y
conda activate mmpbsa
pip install gmx_MMPBSA
gmx_MMPBSA --version
```

### 1.4 Repo and environment record

```bash
git clone <your-remote> brusafe-md && cd brusafe-md
chmod +x scripts/*.sh analysis/*.sh mmpbsa/*.sh
./scripts/00_env.sh
```

This writes `docs/env_<hostname>.txt`. **Commit it from every machine.** When a
judge asks what you ran on, this is the answer, and when two machines disagree
you will see it here first.

### 1.5 Where trajectories live

Default is `runs/` inside the repo. If you have shared storage, point every
machine at it so the duplicate-run guard actually works across boxes:

```bash
echo 'export BRUSAFE_RUNS=/mnt/shared/brusafe/runs' >> ~/.bashrc
```

Without shared storage the guard only protects within one machine, and you must
keep the assignment spreadsheet honest yourself.

---

## 2. Input structures — the manual part

The step-by-step version, including CHARMM-GUI settings and the
`scripts/00_build_inputs.py` helper that assembles the dimers and mutates the
RNA, is in **`docs/TOPOLOGIES.md`**. This section is the reasoning behind it.

Each system expects one of:

```
input/<system>/<system>.pdb        # pdb2gmx route  (apo systems)
input/<system>/gromacs/            # CHARMM-GUI Solution Builder output (RNA)
```

### S1 — WT C/C dimer from 2MS2

The asymmetric unit gives you three quasi-equivalent chains and therefore **two
physically different dimers**: the symmetric C/C pair and the asymmetric A/B
pair. The free coat protein dimer in solution is the C/C-like conformer; A/B is
what assembly and RNA binding induce. Simulating the wrong one silently gives
you an assembly-locked baseline.

Do not take chain IDs on trust. Open the entry, look at the biological assembly,
and identify which pair has **both** FG loops in the extended conformation. Then:

```bash
mkdir -p input/S1_wt_cc_apo
pdb_selchain -C,C 1msc.pdb | pdb_delhetatm | pdb_tidy > input/S1_wt_cc_apo/S1_wt_cc_apo.pdb
# or without pdb-tools:
#   grep '^ATOM' 1msc.pdb | awk 'substr($0,22,1)=="C"' > ...
```

Checklist before you accept a PDB as input:

- [ ] crystallographic waters removed
- [ ] alternate locations resolved to one conformer
- [ ] cryoprotectants / buffer components removed
- [ ] no missing residues in the FG loop (build them if there are)
- [ ] termini decided consciously — `01_prepare.sh` runs `pdb2gmx -ter`
      interactively so you have to choose

### S2 — single-chain dimer

Build with AlphaFold3 or by grafting the linker onto the 2MS2 dimer. **Before
simulating**, superpose the model onto the 2MS2 C/C dimer and record backbone
RMSD:

```bash
gmx confrms -f1 input/S1_wt_cc_apo/S1_wt_cc_apo.pdb \
            -f2 input/S2_sccp_apo/S2_sccp_apo.pdb -o fit.pdb
```

That number goes on the wiki either way. If the two CP subunits sit in a
non-native relative orientation, every S1-vs-S2 comparison downstream measures
an AlphaFold artefact instead of a linker effect, and no amount of sampling
fixes it.

### S4, S8, S9 — RNA-bound systems

Use the CP–operator **cocrystal**. Do not dock an RNA hairpin de novo when a
solved complex exists; reviewers will ask why, and you will not like the answer.
Build your designed pac and the calibration variants by mutating the cocrystal
RNA in place.

For anything containing RNA, prefer CHARMM-GUI Solution Builder over the
`pdb2gmx` route:

```bash
./scripts/01_prepare.sh S4_cp_pacdesign --source charmm-gui
```

Two reasons. `pdb2gmx` residue mapping for RNA under the charmm36 port is
fragile. And CHARMM-GUI handles Mg²⁺ placement properly, which matters more than
it looks — see §7.

---

## 3. S1 end to end

```bash
./scripts/01_prepare.sh     S1_wt_cc_apo --replica 1   # once per replica (own ion placement)
./scripts/02_equilibrate.sh S1_wt_cc_apo 1    # once per replica
```

`01_prepare.sh` prints the index groups it built. You want to see `Backbone`,
`C-alpha`, `SOLU`, `SOLV` in that list — the script now dies if any are missing,
because `04_postprocess.sh` asks for them by name and failing here is much
cheaper than failing three days later.

`02_equilibrate.sh` runs EM → NVT 200 ps restrained → NPT 1 ns restrained →
NPT 2 ns free, then applies the **equilibration gate**: density near 1000 kg/m³,
temperature within 2 K of 310, and density no longer drifting. It refuses to
exit successfully otherwise. A drifting box means the system has not settled
whatever the mean says, and the gate catches that where a mean alone would not.

Take about 20 minutes here. Look at `density.xvg` and `rmsd` yourself. This is
the last cheap moment to notice something is wrong.

---

## 4. Benchmark once, then lock

You need a `prod.tpr` to benchmark against, so make one by hand from the
equilibrated replica:

```bash
cd runs/S1_wt_cc_apo/rep1
gmx grompp -f ../../../mdp/prod.mdp -c npt_free.gro -t npt_free.cpt \
           -p topol.top -n index.ndx -o prod.tpr
cd -

./scripts/bench_gpu.sh runs/S1_wt_cc_apo/rep1/prod.tpr
rm runs/S1_wt_cc_apo/rep1/prod.tpr    # let 03_production regenerate it from the registry
```

That `rm` matters. `03_production.sh` sets `nsteps` from the `ns` column in
`config/systems.tsv`, and the RNA systems run 30 ns rather than 100. A hand-made tpr would
silently override that.

The sweep tries three offload configurations across several thread counts, then
runs two concurrent jobs on one GPU. At ~35k atoms a single run usually does not
saturate a 4070, so **aggregate throughput from two concurrent runs is often
higher than one run alone** — that is the number that sets your schedule.

Then lock it, everywhere:

```bash
echo 'export BRUSAFE_MDRUN_FLAGS="-nb gpu -pme gpu -bonded gpu -update gpu -ntmpi 1 -ntomp 8 -pin on"' >> ~/.bashrc
```

Record the winning string and the date in `docs/DECISIONS.md`. Every replica in
the campaign must use the identical flags or your timings are not comparable.

If one machine benchmarks 25% below the others, check for thermal or power
throttling before blaming the software:

```bash
nvidia-smi -q -d PERFORMANCE | grep -A5 'Clocks Event Reasons'
```

OEM desktops with poor airflow are a real and common source of spread across ten
nominally identical boxes.

---

## 5. Plan across the ten machines

```bash
python3 scripts/_registry.py plan --machines 10 --tier 1 --nsday <your measured number>
```

You get a cost table and a per-machine job list. Tier 1 is about 1,400
ns-equivalent including equilibration; at 400 ns/day per GPU (HMR, 4 fs) that is
well under a day on ten GPUs. Compute is not your constraint;
that is the point of running this early.

Save the job lists and give one to each machine:

```bash
python3 scripts/_registry.py plan --machines 10 --tier 1 --nsday <measured> > docs/PLAN.txt
git add docs/PLAN.txt && git commit -m "campaign plan, tier 1"
```

---

## 6. The production loop

Per replica, on the assigned machine:

```bash
./scripts/01_prepare.sh     S1_wt_cc_apo --replica 2   # skip if build_rep2/ already exists
./scripts/02_equilibrate.sh S1_wt_cc_apo 2
./scripts/03_production.sh  S1_wt_cc_apo 2 23    # 23 = wall-clock hours per chunk
./scripts/04_postprocess.sh S1_wt_cc_apo 2
```

`03_production.sh` is **restartable**. When it hits the wall-clock limit or the
machine reboots, rerun the identical command and it resumes from the last
checkpoint. A simple loop is enough to run a replica to completion overnight:

```bash
until ./scripts/03_production.sh S1_wt_cc_apo 2 23; do sleep 10; done
```

It also writes `PROVENANCE.txt` next to the trajectory — host, GROMACS version,
mdrun flags, seed, start and finish times. That file is what makes a run
defensible six weeks later when nobody remembers which machine did what.

Progress across everything, any time:

```bash
python3 scripts/_registry.py status --runs runs/
```

Run this **daily**. The realistic failure with ten unmanaged desktops is not
running out of compute — it is two boxes running the same replica while a system
never launched at all.

---

## 7. When things fail

**Equilibration gate fails.** Plot `density.xvg`. A slow monotonic drift usually
means the box has not settled: extend `npt_free` and rerun. A sudden jump or a
crash on restraint release usually means a bad starting structure — go back to
§2. Do not force past this gate; production from an unequilibrated box wastes
GPU-days and you will not notice until analysis.

**`already holds a completed run`.** The duplicate guard did its job. Confirm
nobody else is running that replica, then `BRUSAFE_FORCE=1 ./scripts/03_...`
if you genuinely mean to overwrite.

**`gmx dssp` unavailable.** GROMACS is older than 2023. Upgrade, or fall back to
`do_dssp` with an external DSSP binary. Everything else still runs.

**`no 'RNA' group`.** `make_ndx` did not recognise your nucleic acid residue
names. Add the group manually:

```bash
gmx make_ndx -f runs/<sys>/build/solv_ions.gro \
             -n runs/<sys>/build/index.ndx -o runs/<sys>/build/index.ndx
```

**Mg²⁺ in the interface.** `04_postprocess.sh` writes `mg_interface.txt` for
RNA systems (fraction of frames each ion is within 0.5 nm of the RNA, of the
protein, and of both; BOTH > 0 means it sits at the interface), plus the raw
`mindist_mg_rna.xvg` / `mindist_mg_protein.xvg`. Look at it. Mg²⁺ water exchange takes microseconds, so on a 30–100 ns
trajectory the ions sit essentially where you placed them — they do not
equilibrate, and random `genion` placement adds noise rather than realism. An
ion parked in the protein–RNA interface contaminates every M3.5 energy from
that replica. Preferred placement order: crystallographic sites, then
CHARMM-GUI, then `genion` as a documented last resort. This is a good question
to put to Diego and worth writing up for the wiki either way.

---

## 8. Weekly rhythm

- **Daily:** `_registry.py status`. Anything stalled gets restarted the same day.
- **Twice a week:** commit analysis output, `docs/env_*.txt`, `PROVENANCE.txt`.
  Never commit `.xtc` — they are gitignored for a reason. Back trajectories up
  separately; losing them is the one unrecoverable failure here.
- **27 September:** hard stop on new systems. Anything not launched is out of
  scope. Compute was never the constraint; your attention is.

---

## 9. The comparison that carries the claim

A single system's RMSD proves nothing. The claim lives in the comparison:

```bash
./analysis/compare_systems.sh S1_wt_cc_apo S2_sccp_apo
```

This computes the **essential subspace overlap** between the two systems'
principal motions. Over the first 10 eigenvectors, >0.7 supports "the linker
does not alter the dimer's dynamics"; <0.5 means it changed something real. Two
RMSF curves that look similar by eye is not evidence — overlap is.

It also prints inter-replica RMSF agreement within each system. **Read that
first.** A difference between S1 and S2 only means something if it exceeds the
spread among replicas of the same system. Most over-claimed MD results die right
here.

One practical note: S2 has extra linker residues, so the two systems do not have
identical atom counts. Restrict the covariance analysis to CP backbone atoms
present in both before comparing.

---

## 10. M3.5

Only after `04_postprocess.sh` has run on the RNA-bound system. Find your group
numbers, then:

```bash
gmx make_ndx -f runs/S4_cp_pacdesign/build/solv_ions.gro -o /dev/null  # note the numbers, then q

export PROT_GRP=<n> RNA_GRP=<n>
./mmpbsa/run_mmpbsa.sh S4_cp_pacdesign 1 4.0
```

Check the output **in this order**, and the order is not negotiable:

1. `FINAL_DECOMP.dat` — do the per-residue hot spots recover the known
   operator-contact residues of MS2 CP? If not, stop. The setup is wrong and no
   energy from it means anything. This costs a day and gates everything after it.
2. Repeat at `indi` = 1, 2 and 4. Report the spread as part of the error bar.
   A conclusion that flips between 2 and 4 is not a conclusion.
3. Repeat on two disjoint frame windows.
4. Report **ΔΔG between systems only**. Absolute MM/PBSA binding free energies
   are routinely off by 5–20 kcal/mol, worse for charged protein–RNA interfaces.
   Never quote one.
5. ΔΔG = −RT ln(ratio), RT = 0.616 kcal/mol at 310 K — one order of magnitude in
   K_d is about 1.4 kcal/mol.

---

## 11. Before the wiki freeze

- [ ] `docs/env_*.txt` committed from every machine
- [ ] `docs/DECISIONS.md` complete, including the locked mdrun flags and any
      deviation from the defaults
- [ ] `docs/PREDICTIONS.md` timestamped **before** the dPCR readout
- [ ] Systems table: atom counts, box, lengths, replicas, achieved ns/day
- [ ] Block-averaging factor from your own data quoted in the methods section
      (`convergence.py block` prints it) — this preempts the most common
      criticism of student MD analysis
- [ ] An explicit "what we did not validate" section. Writing it yourself is
      always better than having a judge write it for you.
