#!/usr/bin/env bash
# Set up and CHECK this machine for the campaign. Safe to run any number of
# times; on a machine that is already fine it just prints the checklist.
#
#   ./scripts/setup_machine.sh
#
# Prerequisites you do by hand once per PC: Git, a Python 3 with numpy
# (python -m pip install --user numpy), and the repo cloned.
#
# What it does
#   - finds GROMACS in ~/gromacs; if missing, extracts the portable zip from
#     ~/Downloads or from <OneDrive>/BruSafe-MD/software
#   - writes BRUSAFE_GMX into ~/.bashrc (removing any old/broken line)
#   - checks: repo location, Python+numpy, GPU, GROMACS version + CUDA,
#     OneDrive folder; records the machine with 00_env.sh
# and ends with READY, or a list of what to fix.
set -uo pipefail
source "$(dirname "$0")/lib.sh"
set +e   # lib.sh turns errexit on; this script reports problems instead of dying

WANT_GMX="2026.3"
PROBLEMS=()
ok()   { printf '  [ ok ] %s\n' "$*"; }
warn() { printf '  [WARN] %s\n' "$*"; }
bad()  { printf '  [FAIL] %s\n' "$*"; PROBLEMS+=("$*"); }

echo "=== setup check: $(hostname)  $(date '+%Y-%m-%d %H:%M') ==="

# ---- repo location ---------------------------------------------------------
if [ "$REPO" = "$HOME/BruSafe-MD-simulations" ]; then ok "repo at $REPO"
else warn "repo is at $REPO; the plan assumes ~/BruSafe-MD-simulations (mv it there)"; fi
BR="$(git -C "$REPO" rev-parse --abbrev-ref HEAD 2>/dev/null || echo '?')"
[ "$BR" = "claude/gallant-babbage-ydyxr2" ] && ok "branch $BR" \
  || warn "branch is '$BR'; expected claude/gallant-babbage-ydyxr2 (git checkout it)"

# ---- python ----------------------------------------------------------------
if [ -z "$PY" ]; then
  bad "no working Python 3 (install it: winget install -e --id Python.Python.3.12 --scope user)"
elif "$PY" -c 'import numpy' >/dev/null 2>&1; then
  ok "python $("$PY" -c 'import sys; print(sys.version.split()[0])' | tr -d '\r'), numpy $("$PY" -c 'import numpy; print(numpy.__version__)' | tr -d '\r')"
else
  bad "numpy missing: $PY -m pip install --user numpy"
fi

# ---- GPU -------------------------------------------------------------------
if command -v nvidia-smi >/dev/null 2>&1; then
  G="$(nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null | head -1 | tr -d '\r')"
  case "$G" in *4070*) ok "GPU: $G" ;; "") bad "nvidia-smi gave no GPU" ;; *) warn "GPU is not a 4070: $G" ;; esac
else bad "nvidia-smi not found"; fi

# ---- GROMACS ---------------------------------------------------------------
find_gmx() { find "$HOME/gromacs" -name gmx.exe 2>/dev/null | head -1; }
GMX_PATH="$(find_gmx)"
if [ -z "$GMX_PATH" ]; then
  ZIP=""
  for z in "$HOME"/Downloads/gromacs-${WANT_GMX}-win64-cuda*.zip \
           "$HOME"/OneDrive*Monterrey*/BruSafe-MD/software/gromacs-${WANT_GMX}-win64-cuda*.zip \
           "$HOME"/OneDrive*/BruSafe-MD/software/gromacs-${WANT_GMX}-win64-cuda*.zip; do
    [ -f "$z" ] && { ZIP="$z"; break; }
  done
  if [ -z "$ZIP" ]; then
    bad "GROMACS not found: put gromacs-${WANT_GMX}-win64-cuda12.6.3-sm89.zip in ~/Downloads or OneDrive/BruSafe-MD/software and rerun"
  else
    echo "  extracting $ZIP ..."
    mkdir -p "$HOME/gromacs"
    if   [ -x /c/Windows/System32/tar.exe ]; then /c/Windows/System32/tar.exe -xf "$ZIP" -C "$HOME/gromacs"
    elif command -v unzip >/dev/null 2>&1;  then unzip -q -o "$ZIP" -d "$HOME/gromacs"
    else bad "no tool to extract the zip (tar.exe / unzip)"; fi
    GMX_PATH="$(find_gmx)"
  fi
fi
if [ -n "$GMX_PATH" ]; then
  # one clean line in ~/.bashrc, whatever was there before
  touch "$HOME/.bashrc"
  sed -i '/^export BRUSAFE_GMX=/d' "$HOME/.bashrc"
  echo "export BRUSAFE_GMX=\"$GMX_PATH\"" >> "$HOME/.bashrc"
  V="$("$GMX_PATH" --version 2>/dev/null | tr -d '\r')"
  VER="$(printf '%s\n' "$V" | awk -F': *' '/^GROMACS version/{print $2; exit}')"
  GPUS="$(printf '%s\n' "$V" | awk -F': *' '/^GPU support/{print $2; exit}')"
  if [ "$VER" = "$WANT_GMX" ] && [ "$GPUS" = "CUDA" ]; then ok "GROMACS $VER ($GPUS) at $GMX_PATH"
  else bad "GROMACS is '$VER' / GPU '$GPUS' at $GMX_PATH; the campaign needs $WANT_GMX with CUDA"; fi
fi

# ---- OneDrive --------------------------------------------------------------
OD=""
for od in "$HOME"/OneDrive*Monterrey* "$HOME"/OneDrive*; do [ -d "$od" ] && { OD="$od"; break; }; done
if [ -n "$OD" ]; then ok "OneDrive: $OD"
else warn "no OneDrive folder: results must be uploaded from the browser (pack/<host>)"; fi

# ---- record this machine -----------------------------------------------------
if [ "${#PROBLEMS[@]}" -eq 0 ]; then
  BRUSAFE_GMX="$GMX_PATH" "$REPO/scripts/00_env.sh" > /dev/null 2>&1 \
    && ok "recorded docs/env_$(hostname).txt" || warn "00_env.sh failed; run it by hand"
fi

echo
if [ "${#PROBLEMS[@]}" -eq 0 ]; then
  echo "READY: $(hostname). Open a NEW terminal (so ~/.bashrc is read) before the first run."
else
  echo "NOT READY: ${#PROBLEMS[@]} problem(s):"; printf '  - %s\n' "${PROBLEMS[@]}"
  exit 1
fi
