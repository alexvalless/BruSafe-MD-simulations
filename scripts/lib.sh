
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MDP="$REPO/mdp"
CFG="$REPO/config"
RUNS="${BRUSAFE_RUNS:-$REPO/runs}"

# Locked mdrun configuration. Set by bench_gpu.sh
# Every replica in the campaign must use the identical string.
MDRUN_FLAGS="${BRUSAFE_MDRUN_FLAGS:--nb gpu -pme gpu -bonded gpu -update gpu -ntmpi 1 -ntomp 8 -pin on}"

# Solvent/ion residue names excluded from the SOLU group.
# Covers the MacKerell GROMACS port (SOL/POT/CLA/MG) and CHARMM-GUI (TIP3).
SOLVENT_RESNAMES="SOL TIP3 TIP3P POT CLA SOD MG MGA ZN2 CAL"

log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die()  { printf '[FATAL] %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "missing required command: $1"; }

# Read one field from config/systems.tsv.
# usage: sysfield <system_name> <column_name>
sysfield() {
  python3 "$REPO/scripts/_registry.py" field "$1" "$2"
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
