# Running on Windows without administrator rights (Git Bash)

The scripts are plain bash and run unchanged on Linux and WSL. This page covers laboratory machines
on which WSL cannot be enabled. Everything below installs into the user's own folder; no
administrator rights are needed.

If WSL **is** available (`wsl --status` answers without error), it is preferred, and
`docs/TUTORIAL.md` applies as written.

## Requirements for each machine

| Component | Installation without administrator rights |
| --- | --- |
| NVIDIA driver | Already installed on the RTX 4070 machines. Check with `nvidia-smi` |
| bash + coreutils | **PortableGit** from git-scm.com (self-extracting folder; run `git-bash.exe`) |
| Python 3 + numpy | python.org installer, option *Install for me only*; then `python -m pip install --user numpy` |
| GROMACS (CUDA) | Portable `gmx.exe` folder, built once on a machine with administrator rights (see below) |
| CHARMM36m port | Unzip `charmm36-jul2022.ff` into a folder, for example `C:\Users\<user>\ff` |

**One binary for the whole campaign.** GROMACS is built once and the same folder is copied to every
machine. Different GROMACS versions across replicas invalidate the error bars (see `TUTORIAL.md`,
section 1.1).

### Building the portable GROMACS (once, on a machine with administrator rights)

Visual Studio Build Tools (C++ workload), CMake and the CUDA Toolkit are installed, and the build is
run from an *x64 Native Tools* prompt:

```
cmake -S gromacs-2024.x -B build -G Ninja -DCMAKE_BUILD_TYPE=Release ^
      -DGMX_GPU=CUDA -DGMX_FFT_LIBRARY=fftpack -DGMX_HWLOC=OFF ^
      -DCMAKE_INSTALL_PREFIX=C:\gromacs
cmake --build build --target install
```

The CUDA runtime DLLs on which `gmx.exe` depends (`cudart64_*.dll`, `cufft64_*.dll`, from the `bin`
folder of the CUDA Toolkit) are copied into `C:\gromacs\bin`. The folder is then portable and can be
copied to each laboratory machine.

## One-time setup per machine (inside Git Bash)

```bash
# 1. line endings: prevent git from converting the scripts to CRLF
git config --global core.autocrlf false

# 2. clone (the repository's .gitattributes also forces LF)
git clone https://github.com/alexvalless/BruSafe-MD-simulations.git
cd BruSafe-MD-simulations

# 3. environment: adjust the paths, then reopen Git Bash
cat >> ~/.bashrc << 'RC'
export BRUSAFE_GMX="/c/gromacs/bin/gmx.exe"   # or wherever the folder is located
export GMXLIB="$HOME/ff"                      # folder containing charmm36-jul2022.ff
# export BRUSAFE_PYTHON="/c/Users/<user>/AppData/Local/Programs/Python/Python312/python.exe"
RC

# 4. record the state of the machine and check that every line is as expected
./scripts/00_env.sh
```

`00_env.sh` should show `os : windows`, a CUDA-enabled GROMACS, the RTX 4070, a Python 3 version,
numpy, and the charmm36 folder under `GMXLIB`.

## What the scripts handle automatically on Windows

`scripts/lib.sh` detects Git Bash and:

- routes every `gmx` call to `BRUSAFE_GMX` if it is set;
- converts `GMXLIB` to a `C:/...` path that native Windows programs understand;
- finds a working Python 3 and skips the Microsoft Store `python3` stub;
- removes the `\r` that Windows Python adds to its output (otherwise `yes` becomes `yes\r` and
  registry lookups fail silently);
- omits `-pin on`, which GROMACS does not support on Windows (a speed effect only; the results are
  unaffected);
- runs `pdb2gmx` through `winpty`, so that its interactive questions about termini work in the
  Git Bash window.

## Running

The commands are the same as on other platforms. Once per system, before the first night, the
structures are prepared (the RNA systems come from CHARMM-GUI and do not need `GMXLIB`; S1, S2 and S6
go through `pdb2gmx` and do):

```bash
./scripts/01_prepare.sh S4_cp_pacdesign --source charmm-gui --replica 2
./scripts/01_prepare.sh S8_pac_free --source charmm-gui --replica 1
```

Ion placement is seeded for each replica (`build_rep<N>`), so a given replica builds identically on
every machine with the same GROMACS build, while different replicas have different ion positions as
well as different velocities. `run_slot.sh` performs this step automatically; to do it by hand,
`--replica N` is used. For overnight runs:

```bash
./scripts/run_night.sh S1_wt_cc_apo 1 11
```

- **The screen should be locked (Win+L), not signed out.** Signing out terminates the run. If it
  terminates anyway, the same command resumes from the checkpoint.
- The machine should not be set to sleep or to restart for updates overnight.
- In the morning the replica folder from `runs/` is **copied** (not moved) to the shared drive.
  `prod.cpt` must remain on the machine for the next night.

## Known limitations

- **MM/PBSA (M3.5) does not run on native Windows.** gmx_MMPBSA depends on AmberTools.
  `mmpbsa/run_mmpbsa.sh` is run on a Linux, WSL or macOS machine, using the trajectories produced
  here.
- If `BRUSAFE_MDRUN_FLAGS` is overridden on Windows, `-pin on` should be left out.
