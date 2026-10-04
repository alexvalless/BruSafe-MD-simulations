# Building the starting structures and topologies

Everything in this file is done **once per system, on one machine**, and the
result is committed. The school machines then only run equilibration and
production (`02_equilibrate.sh`, `03_production.sh`), so they never repeat this
work and every replica of a system starts from the identical topology.

```
input/<system>/<system>.pdb        built by scripts/00_build_inputs.py      (commit)
input/<system>/gromacs/            CHARMM-GUI output, RNA systems only      (commit)
runs/<system>/build/               01_prepare.sh: solvated, ionised, HMR,   (copy to
                                   index.ndx                                 the GPUs)
```

---

## 0. Before you start

- Same GROMACS version on the build machine and on every GPU machine
  (`gmx --version`). A `.tpr` is not portable between versions, and the
  scripts run `grompp` on the GPU machine, so the version must match there.
- CHARMM36m GROMACS port from the MacKerell lab unpacked under `$GMXLIB`, as
  `charmm36-jul2022.ff` (see `docs/TUTORIAL.md` §1.2). Another directory
  name: `export BRUSAFE_FF=<name>` and add a line to `docs/DECISIONS.md`.
- Python 3.9+, no extra packages, for `scripts/00_build_inputs.py` and
  `scripts/hmr_top.py`.

---

## 1. Fill in the registry

`config/systems.tsv`, column `source_pdb`:

| value | meaning |
|---|---|
| `1MSC` (any PDB ID) | downloaded to `input/_pdb/1MSC.pdb` on first use |
| `model` | your own model at `input/<system>/model.pdb` (AlphaFold, linker graft) |
| `from:<system>` | take chains from another system's built PDB |
| `cocrystal` | placeholder for a PDB ID (the registry now uses `1ZDH`) |

For the RNA-bound systems use **one** MS2 coat protein–operator structure for
S3, S4, S5 and S7, so that the only thing that changes between them is the RNA
sequence. No docking is needed: MS2 capsid crystals soaked with the 19-nt
operator hairpin put one hairpin on the A/B dimer of the asymmetric unit
(Valegård et al. 1994, 1997; Grahn et al.). The registry uses **1ZDH**
(wild-type operator, 2.7 Å). Alternatives worth knowing:

| entry | what it is | use |
|---|---|---|
| 1ZDH | WT operator, 2.7 Å | starting structure for S3–S5, S7 |
| 1ZDI | operator variant, 2.7 Å (reported as the C(−5) variant — confirm) | check S7_cal_u5c: mutate 1ZDH in silico, compare with this crystal |
| 2BU1 | 5-bromo-U at −5, 2.2 Å | higher resolution; `build` reverts 5BU to U |

Confirm chain letters and resolved nucleotides with `show` before building,
and record the choice in `docs/DECISIONS.md`.

`config/rna_designs.tsv` holds the RNA sequence of each RNA-bound system.
`native` keeps the crystal sequence; `TODO` blocks the build on purpose.

---

## 2. Look before you build

```bash
python3 scripts/00_build_inputs.py show 1MSC
```

It prints, per chain: sequence, residue range, chain breaks and missing heavy
atoms, then the **symmetry partners** of each chain from the BIOMT operators,
ranked by number of contacts (illustrative output):

```
C <- C op7   312 contacts  rotation 180.0 deg      <- the C/C dimer (two-fold)
A <- B op1   298 contacts  rotation   0.0 deg      <- A/B dimer, both in the file
```

What to check here, because nothing downstream can:

- **S1 must be the C/C pair, S6 the A/B pair.** The C/C partner has to be a
  180° operator. The FG loops (residues ~66–82) should both be extended in C/C;
  in A/B one is folded back. Open the result in PyMOL/VMD and look.
- **No `BREAK` in the FG loop or anywhere else.** A break means missing residues;
  model them (MODELLER, or CHARMM-GUI / pdbfixer for short loops) before going
  on. pdb2gmx would otherwise join the ends with a bond across the gap.
- For the cocrystal: which chain letter is the RNA, how many nucleotides are
  resolved, and whether there are crystallographic Mg²⁺ next to it.

---

## 3. Build the PDB

```bash
python3 scripts/00_build_inputs.py build S1_wt_cc_apo
python3 scripts/00_build_inputs.py build S6_wt_ab_apo
python3 scripts/00_build_inputs.py build S3_cp_operator
python3 scripts/00_build_inputs.py build S4_cp_pacdesign      # after filling rna_designs.tsv
...
```

What `build` does: first model only, one alternate location per residue,
MSE→MET, hydrogens/waters/other HETATM removed (crystallographic Mg²⁺ within 6 Å
kept for RNA systems, `--no-mg` to drop them), second and later chains placed with
the symmetry operator that gives the most contacts, RNA mutated in place. Output
chains are renamed `A, B, …` (protein), `R, S, …` (RNA), `M` (Mg²⁺). Every
choice it made is written as `REMARK 1` lines at the top of the PDB, so the file
documents itself.

The RNA mutation keeps the sugar-phosphate backbone and the base atoms shared by
the old and new base. For purine↔pyrimidine swaps the glycosidic atoms are
renamed (N9→N1, C4→C2, C8→C6) so the new base is rebuilt in the same *anti*
orientation. The atoms it removes are listed as "missing heavy atoms", and
**CHARMM-GUI rebuilds them** in the next step (pdb2gmx cannot, which is one more
reason the RNA systems go through CHARMM-GUI).

---

## 4a. Apo systems (S1, S2, S6): pdb2gmx

```bash
./scripts/01_prepare.sh S1_wt_cc_apo
```

pdb2gmx asks for the termini of each chain (`-ter`): choose **NH3+** and
**COO-** for the native CP termini (Ala1 / Tyr129). Histidine tautomers are
assigned by pdb2gmx from the hydrogen-bond network; if a His sits at the
RNA-binding face, check the choice in `runs/<system>/build/topol*.itp`.

Then the script solvates (rhombic dodecahedron, 1.0 nm), adds 150 mM KCl,
**applies HMR to the topology** and builds `index.ndx`. Expected output includes

```
topol_Protein_chain_A.itp:Protein_chain_A: 1xxx hydrogens repartitioned (total mass ... unchanged)
HMR: ...
```

## 4b. RNA systems (S3–S5, S7, S8, S9): CHARMM-GUI Solution Builder

At charmm-gui.org → Input Generator → **Solution Builder**:

1. **Upload** `input/<system>/<system>.pdb` (PDB format). Select all chains:
   the two protein segments, the RNA, and the Mg²⁺ (HETATM) if present.
2. **PDB manipulation.** Confirm the nucleic acid chain is read as **RNA**.
   RNA termini: 5′-OH / 3′-OH (no terminal phosphate) unless your construct
   needs one. Protein termini: NTER / CTER. No mutations here, they were done in
   step 3. Missing atoms (the mutated bases) are built automatically.
3. **Water box.** Octahedral, *fit to protein size*, edge distance **10 Å**.
   (S8, the free hairpin: 12 Å.)
4. **Ions.** KCl **0.15 M**, neutralising. Mg²⁺: if crystallographic Mg²⁺ were
   kept, do not add more unless `mg_count` in the registry asks for it; if not,
   add MgCl₂ so the count matches `mg_count`, and write the placement method in
   `docs/DECISIONS.md` (see the Mg²⁺ entry there).
5. **Input generation.** Force field **CHARMM36m**, output **GROMACS**,
   temperature **310 K**. Leave **hydrogen mass repartitioning unticked**:
   `hmr_top.py` applies it in the next step exactly as for the apo systems. (If
   you do tick it, `hmr_top.py` detects heavy hydrogens and skips itself; never
   both.)
6. Download, unzip, and copy the `gromacs/` folder to `input/<system>/gromacs/`
   (it must contain `step3_input.gro`, `topol.top`, `toppar/`).

```bash
./scripts/01_prepare.sh S3_cp_operator --source charmm-gui
```

If it warns that there is no `RNA` group, add one with `gmx make_ndx` as it
tells you; `04_postprocess.sh` and MM/PBSA need it.

### Single-chain CP (S2, S9) and the free hairpin (S8)

- **S2 / S9:** save the model as `input/<system>/model.pdb`, run `build`, then
  superpose on S1 and record the backbone RMSD before simulating
  (`docs/TUTORIAL.md` §2). S9 needs the RNA placed as in the cocrystal: build
  S4 first and superpose the scCP model onto its protein chains.
- **S8:** either a ViennaRNA-guided 3D model as `model.pdb`, or the bound
  conformation from S4 (`source_pdb = from:S4_cp_pacdesign`, `chains = R`), which
  is the usual starting point for a pre-organisation penalty. Record the choice.
  Build with `BRUSAFE_BOX_PAD=1.2` / 12 Å.

---

## 5. Check the build, then hand it to the GPU machines

On the build machine, before anything goes to a GPU:

```bash
cd runs/S1_wt_cc_apo/build
head -1 topol*.itp toppar/*.itp 2>/dev/null | grep HMR    # marker present
gmx grompp -f ../../../mdp/prod.mdp -c solv_ions.gro -p topol.top -n index.ndx -o /tmp/check.tpr
```

`grompp` must finish with no WARNING. One NOTE is expected with HMR at 4 fs:
*"bond … has an estimated oscillational period … less than 10 times the time
step"*. Anything else, stop and look.

Commit `input/<system>/` and copy `runs/<system>/build/` to the GPU machines
(or point `BRUSAFE_RUNS` at shared storage). There, per replica:

```bash
./scripts/02_equilibrate.sh S1_wt_cc_apo 1
until ./scripts/03_production.sh S1_wt_cc_apo 1 23; do sleep 10; done
```

---

## Atom counts and expected cost

Record these in the systems table on the wiki. Rough expectations, rhombic
dodecahedron / octahedron at 1.0 nm:

| system | atoms (approx.) |
|---|---|
| CP dimer apo (S1, S6) | 40–60 k |
| CP dimer + 19-nt hairpin (S3–S5, S7) | 60–100 k (the hairpin protrudes) |
| free hairpin, 1.2 nm (S8) | 15–25 k |

Measure ns/day once with `bench_gpu.sh` and plan with
`python3 scripts/_registry.py plan --machines <n> --tier 1 --nsday <measured>`.
