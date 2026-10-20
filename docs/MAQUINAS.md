# Laboratory machines: step-by-step guide

All commands are typed in the **MSYS2 UCRT64** terminal (Start menu → "MSYS2 UCRT64"), **one at a
time**, exactly as written, without changes (except the slot number in step 8). The first time on a
PC takes about 25 minutes. After that, the job runs unattended (1 to 2 hours per job).

The following is already settled and does not need to be changed: the GROMACS version is **2026.3 with
CUDA** (the same on every machine), the OneDrive path is stored in the scripts, and the results are
copied automatically when a job finishes.

---

## 1. GPU and Git

```bash
git --version; nvidia-smi --query-gpu=name --format=csv,noheader
```

The output should be the Git version and `NVIDIA GeForce RTX 4070`. If the GPU is not listed, the PC
cannot be used.

## 2. Install Python (without administrator rights)

```bash
winget.exe install -e --id Python.Python.3.12 --scope user --accept-package-agreements --accept-source-agreements
```

If `winget` does not exist or fails, use the installer from python.org (3.12, 64-bit) with **"Install
for me only"** and "Add python.exe to PATH" ticked.

**Close the terminal and open it again.** Then:

```bash
python --version
```

The output should be `Python 3.12.x`. Use `python` and not `python3`, which is a Microsoft Store
shortcut that does not work (the project scripts already ignore it).

## 3. Install numpy

```bash
python -m pip install --user numpy
```
```bash
python -c "import numpy; print(numpy.__version__)"
```

The output should be a version number (for example `2.5.3`).

## 4. Clone the repository

```bash
cd ~ && git clone https://github.com/alexvalless/BruSafe-MD-simulations.git && cd BruSafe-MD-simulations
```

## 5. The GROMACS zip

```bash
ls -lh "/c/Users/A01563079/OneDrive - Instituto Tecnologico y de Estudios Superiores de Monterrey/BruSafe-MD/software/"
```

The file `gromacs-2026.3-win64-cuda12.6.3-sm89.zip` of **about 242 MB** should appear.

- **It appears but shows no size, or 0:** OneDrive has only the file name and not the file. In
  Windows Explorer, right-click the zip → **"Always keep on this device"** and wait for the solid
  green mark.
- **The folder does not exist:** OneDrive is not running on that PC. Download the zip to
  `Downloads` while signed in to GitHub (and sign out again afterwards):
  https://github.com/alexvalless/BruSafe-MD-simulations/actions/runs/37343539360/artifacts/11363541517
- Integrity check (optional): the SHA-256 sum should be
  `ed5fd8e86a8733433a47cba8f816955070d20369ef10fbff6e77a6f38760a940`.

## 6. Prepare and verify the machine

```bash
./scripts/setup_machine.sh
```

The script extracts GROMACS, sets up `~/.bashrc` and checks everything. It should end with
`READY: CHI-12104-PC-XX`. If it prints `NOT READY`, it lists what is missing. **Close the terminal
and open a new one** before continuing.

## 7. Check the plan

The jobs are listed in `docs/PLAN.txt`, one line per job: slot (`gpu<N>`), system, replica, length in
ns and the purpose of the job. Each PC takes one slot number and runs all the jobs of that slot. Write
down which PC took which number so that no work is repeated.

The command sheets `docs/Comandos_S2_S9.docx` and `docs/Comandos_S3.docx` give the launch commands for
the current plan; this guide describes the general procedure.

## 8. Launch (Terminal 1)

With the slot number in place of `N` (and `16` as the time budget in hours):

```bash
cd ~/BruSafe-MD-simulations && git pull --no-rebase && ./scripts/run_slot.sh N 16 > slot_N.log 2>&1
```

This terminal stays **busy and silent** until the job finishes.

## 9. Follow the progress (Terminal 2, a separate UCRT64 window)

The same command works on every machine and refreshes every 30 seconds:

```bash
cd ~/BruSafe-MD-simulations && ./scripts/progress.sh watch
```

| Stage | Duration | What is shown |
|---|---|---|
| build + equilibration | ≈10 min | stage `equil`; in `slot_N.log`: `building`, `equilibration` |
| production | depends on the system (about 1 h for S8, 1.3 h for S4) | `prod`, ns at the last checkpoint (every 15 min), size of `prod.xtc` and how long ago it grew |
| post-processing + packing | a few minutes | `post-processing`, `results packed` |

Leaving Terminal 2 with **Ctrl+C** is safe: it does not affect the job.

## 10. Rules while a job is running

1. **Never press Ctrl+C in Terminal 1.** (This is what stopped PC-02 the first time.)
2. Do not close that window.
3. Lock the screen with **Win+L**; do not sign out and do not let the machine sleep.
4. After a power cut or a restart, open UCRT64 again and **repeat the same command of step 8**. The job resumes from the last checkpoint and does not repeat finished work.

## 11. When it finishes

The prompt of Terminal 1 returns. Check with:

```bash
tail -3 ~/BruSafe-MD-simulations/slot_N.log
```

(with the slot number in place of `N`). It should say `all jobs of slot N completed`. The results are
left in the OneDrive folder `BruSafe-MD/<name-of-the-PC>`. Wait until the files have the **green mark**
before signing out, and sign out of OneDrive and GitHub in the browser (the computers are shared).

---

## If something fails

| Symptom | What to do |
|---|---|
| `missing required command: gmx` | `git pull --no-rebase` (it already fixes this) and relaunch; or `source ~/.bashrc` |
| `Received the INT signal` in the log | Someone pressed Ctrl+C. Relaunch the same command of step 8 |
| `NOT READY` | Read the `[FAIL]` line, correct it and repeat `setup_machine.sh` |
| `python3: Python was not found` | Use `python` (step 2) |
| The zip shows no size | Step 5: "Always keep on this device" |
| The trajectory does not grow for more than 10 min | Report it with the slot number |
| `NO CHARMM-GUI build for ...` | That system is not ready yet; it is not an error of the PC |

To ask for help, copy the name of the PC, the slot number and the output of

```bash
grep -E "^\[|Fatal|FATAL|rror" ~/BruSafe-MD-simulations/slot_N.log | tail -15
```
