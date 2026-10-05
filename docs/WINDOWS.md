# Running on Windows without admin rights (Git Bash)

The scripts are plain bash and run unchanged on Linux and WSL. This page covers
lab machines where WSL cannot be enabled. Everything below installs into your
own user folder; nothing needs administrator rights.

If WSL **is** available (`wsl --status` answers without error), prefer it and
follow `docs/TUTORIAL.md` as written.

## What each machine needs

| Piece | How, without admin |
| --- | --- |
| NVIDIA driver | Already installed on the RTX 4070 machines. Check: `nvidia-smi` |
| bash + coreutils | **PortableGit** from git-scm.com (self-extracting folder, run `git-bash.exe`) |
| Python 3 + numpy | python.org installer, choose *Install for me only*; then `python -m pip install --user numpy` |
| GROMACS (CUDA) | Portable `gmx.exe` folder, built once on a machine that has admin (see below) |
| CHARMM36m port | Unzip `charmm36-jul2022.ff` into a folder, e.g. `C:\Users\<you>\ff` |

**One binary for the whole campaign.** Build GROMACS once and copy the same
folder to every machine. Mixed GROMACS versions across replicas invalidate the
error bars (see `TUTORIAL.md` 1.1).

### Building the portable GROMACS (once, on a machine with admin)

Install Visual Studio Build Tools (C++ workload), CMake and the CUDA Toolkit,
then from an *x64 Native Tools* prompt:

```
cmake -S gromacs-2024.x -B build -G Ninja -DCMAKE_BUILD_TYPE=Release ^
      -DGMX_GPU=CUDA -DGMX_FFT_LIBRARY=fftpack -DGMX_HWLOC=OFF ^
      -DCMAKE_INSTALL_PREFIX=C:\gromacs
cmake --build build --target install
```

Copy the CUDA runtime DLLs that `gmx.exe` depends on (`cudart64_*.dll`,
`cufft64_*.dll`, from the CUDA Toolkit `bin` folder) into `C:\gromacs\bin`.
That folder is now portable: copy it to each lab machine.

## One-time setup per machine (inside Git Bash)

```bash
# 1. line endings: never let git turn the scripts into CRLF
git config --global core.autocrlf false

# 2. clone (the repo's .gitattributes also forces LF)
git clone https://github.com/alexvalless/BruSafe-MD-simulations.git
cd BruSafe-MD-simulations

# 3. environment -- adjust paths, then reopen Git Bash
cat >> ~/.bashrc << 'RC'
export BRUSAFE_GMX="/c/gromacs/bin/gmx.exe"   # or wherever the folder lives
export GMXLIB="$HOME/ff"                      # folder containing charmm36-jul2022.ff
# export BRUSAFE_PYTHON="/c/Users/<you>/AppData/Local/Programs/Python/Python312/python.exe"
RC

# 4. record the machine state and check every line says what you expect
./scripts/00_env.sh
```

`00_env.sh` must show: `os : windows`, a CUDA-enabled GROMACS, the RTX 4070,
a Python 3 version, numpy, and the charmm36 folder under `GMXLIB`.

## What the scripts handle for you on Windows

`scripts/lib.sh` detects Git Bash and:

- routes every `gmx` call to `BRUSAFE_GMX` if set;
- converts `GMXLIB` to a `C:/...` path that native Windows programs understand;
- finds a real Python 3, skipping the Microsoft Store `python3` stub;
- strips the `\r` that Windows Python adds to its output (otherwise `yes`
  becomes `yes\r` and registry lookups silently fail);
- drops `-pin on`, which GROMACS does not support on Windows (speed only,
  results are unaffected);
- runs `pdb2gmx` through `winpty` so its interactive termini questions work
  in the Git Bash window.

## Running

Same commands as everywhere else. For overnight work:

```bash
./scripts/run_night.sh S1_wt_cc_apo 1 11
```

- **Lock the screen (Win+L), do not sign out.** Signing out kills the run.
  If it dies anyway, rerun the same command: it resumes from the checkpoint.
- Check that the machine is not set to sleep or reboot for updates overnight.
- In the morning, **copy** (do not move) the replica folder from `runs/` to
  the shared drive. `prod.cpt` must stay on the machine for the next night.

## Known limits

- **MM/PBSA (M3.5) does not run on native Windows.** gmx_MMPBSA depends on
  AmberTools. Run `mmpbsa/run_mmpbsa.sh` on a Linux, WSL or macOS machine
  using the trajectories produced here.
- If you override `BRUSAFE_MDRUN_FLAGS` on Windows, leave out `-pin on`.
