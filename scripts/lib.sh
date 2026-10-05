#!/usr/bin/env bash
# Shared helpers. Source this, do not execute it.
# Runs on Linux, WSL, macOS and Windows Git Bash (see docs/WINDOWS.md).

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MDP="$REPO/mdp"
CFG="$REPO/config"
RUNS="${BRUSAFE_RUNS:-$REPO/runs}"

log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die()  { printf '[FATAL] %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "missing required command: $1"; }

# ---- platform ---------------------------------------------------------------
case "$(uname -s)" in
  MINGW*|MSYS*|CYGWIN*) BRUSAFE_OS=windows ;;
  Darwin)               BRUSAFE_OS=mac ;;
  *)                    BRUSAFE_OS=linux ;;
esac

# Optional explicit gmx binary, e.g. a portable Windows build on a USB stick:
#   export BRUSAFE_GMX="/c/gromacs/bin/gmx.exe"
# Every `gmx ...` call in the scripts then goes to that binary.
if [ -n "${BRUSAFE_GMX:-}" ]; then
  [ -x "$BRUSAFE_GMX" ] || die "BRUSAFE_GMX=$BRUSAFE_GMX is not an executable"
  gmx() { "$BRUSAFE_GMX" "$@"; }
fi
GMX_BIN="${BRUSAFE_GMX:-gmx}"

# Native Windows programs need GMXLIB as a Windows path (C:/...), not /c/...
if [ "$BRUSAFE_OS" = windows ] && [ -n "${GMXLIB:-}" ] && command -v cygpath >/dev/null 2>&1; then
  GMXLIB="$(cygpath -m "$GMXLIB")"; export GMXLIB
fi

# ---- python -----------------------------------------------------------------
# Picks the first interpreter that is real Python 3. On Windows `python3` is
# often the Microsoft Store stub, which exists but fails -- it is skipped.
# Override with BRUSAFE_PYTHON=/path/to/python.
find_python() {
  local c
  for c in "${BRUSAFE_PYTHON:-}" python3 python py; do
    [ -n "$c" ] || continue
    command -v "$c" >/dev/null 2>&1 || continue
    "$c" -c 'import sys; sys.exit(0 if sys.version_info[0] == 3 else 1)' \
      >/dev/null 2>&1 && { printf '%s\n' "$c"; return 0; }
  done
  return 1
}
PY="$(find_python)" || PY=""

# Run Python with CR stripped from its output. Python on Windows prints \r\n,
# and a stray \r turns "yes" into "yes\r" and breaks every comparison.
pyrun() {
  [ -n "$PY" ] || die "no working Python 3 found (set BRUSAFE_PYTHON)"
  "$PY" "$@" | tr -d '\r'
}

# ---- portable coreutils -----------------------------------------------------
ncpu()    { nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 8; }
now_iso() { date -Is 2>/dev/null || date +%Y-%m-%dT%H:%M:%S%z; }
abspath() { (cd "$(dirname "$1")" && printf '%s/%s\n' "$(pwd)" "$(basename "$1")"); }

# ---- mdrun ------------------------------------------------------------------
# Thread pinning is not supported by GROMACS on Windows; it only affects speed,
# never results, so it is dropped there (docs/DECISIONS.md, 2026-10-03).
if [ "$BRUSAFE_OS" = windows ]; then PIN_FLAGS=""; else PIN_FLAGS="-pin on"; fi

# Locked mdrun configuration. Set by bench_gpu.sh
# Every replica in the campaign must use the identical string.
MDRUN_FLAGS="${BRUSAFE_MDRUN_FLAGS:--nb gpu -pme gpu -bonded gpu -update gpu -ntmpi 1 -ntomp 8 $PIN_FLAGS}"

# Solvent/ion residue names excluded from the SOLU group.
# Covers the MacKerell GROMACS port (SOL/POT/CLA/MG) and CHARMM-GUI (TIP3).
SOLVENT_RESNAMES="SOL TIP3 TIP3P POT CLA SOD MG MGA ZN2 CAL"

# RNA residue names: CHARMM-GUI writes ADE/CYT/GUA/URA, which make_ndx does
# not recognise as RNA, so 01_prepare.sh builds the RNA groups itself.
RNA_RESNAMES="ADE CYT GUA URA A C G U RA RC RG RU"

# Number of a named group in an index file (0-based, as gmx tools count).
# usage: ndxgroup <index.ndx> <name>   -- prints nothing if absent
ndxgroup() {
  awk -v g="$2" '/^\[/ {gsub(/[][ \r]/, ""); if ($0 == g) {print n; exit}; n++}' "$1"
}

# Read one parameter from an mdp file (first match, comments stripped).
# usage: mdpval <file.mdp> <key>
mdpval() {
  awk -v k="$2" -F'=' '{sub(/\r$/, ""); sub(/;.*/, ""); gsub(/[ \t]/, "", $1)}
       $1 == k {gsub(/[ \t]/, "", $2); print $2; exit}' "$1"
}

# grompp for the 4 fs (HMR) stages. HMR leaves bonds without hydrogens alone,
# so C=O (period ~19 fs) trips grompp's "period < 5 x dt" WARNING at 4 fs even
# though HMR at 4 fs is standard practice (Hopkins et al., JCTC 2015). That one
# warning is accepted; ANY other warning still stops the run.
# usage: grompp_hmr <normal gmx grompp arguments>
grompp_hmr() {
  local out n_all n_ok
  if ! out="$(gmx grompp "$@" -maxwarn 10 2>&1)"; then
    printf '%s\n' "$out" >&2; die "grompp failed"
  fi
  printf '%s\n' "$out" >&2
  n_all="$(grep -c '^WARNING [0-9]' <<<"$out" || true)"
  n_ok="$(grep -A3 '^WARNING [0-9]' <<<"$out" | grep -c 'estimated oscillational period' || true)"
  [ "$n_all" -eq "$n_ok" ] || die "grompp gave warnings other than the expected HMR bond-period one -- read them above"
  [ "$n_all" -eq 0 ] || log "accepted $n_all HMR bond-period warning(s) (dt = 4 fs, see docs/DECISIONS.md)"
}

# Read one field from config/systems.tsv.
# usage: sysfield <system_name> <column_name>
sysfield() {
  pyrun "$REPO/scripts/_registry.py" field "$1" "$2"
}

# Deterministic per-replica seed: the same system+replica always yields the
# same seed, so any equilibration can be reproduced from the repo alone.
replica_seed() {
  printf '%s_rep%s' "$1" "$2" | cksum | awk '{print $1 % 2000000000}'
}

# Refuse to clobber a finished production run unless BRUSAFE_FORCE=1.
# Prevents the classic "two machines ran the same replica" failure.
guard_existing() {
  local d="$1"
  if [ -f "$d/prod.gro" ] && [ "${BRUSAFE_FORCE:-0}" != "1" ]; then
    die "$d already holds a completed run (set BRUSAFE_FORCE=1 to override)"
  fi
}

# The pdb2gmx route needs the CHARMM36m port under GMXLIB. Fail early with a
# readable message instead of a GROMACS "force field not found" deep in a run.
check_ff() {
  local ff="${1:-charmm36-jul2022}"
  [ -n "${GMXLIB:-}" ] || die "GMXLIB is not set -- see docs/TUTORIAL.md section 1.2"
  [ -d "$GMXLIB/$ff.ff" ] || die "force field $ff.ff not found in GMXLIB=$GMXLIB"
}
