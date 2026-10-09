#!/usr/bin/env bash
# One unattended night for ONE replica: equilibrate if needed, then run
# production until the time budget is spent. Rerun the identical command the
# next night and it resumes from the last checkpoint.
#
#   ./scripts/run_night.sh <system_name> <replica_number> [total_hours]
#
# Example (leave at 19:00, back at 07:00 -> 11 h budget):
#   ./scripts/run_night.sh S1_wt_cc_apo 1 11
#
# Requires 01_prepare.sh to have run for this system. Do that during the day:
# pdb2gmx asks about termini interactively, so it cannot run unattended.
set -euo pipefail
source "$(dirname "$0")/lib.sh"
need gmx

SYS="${1:?usage: run_night.sh <system> <replica> [hours]}"
REP="${2:?usage: run_night.sh <system> <replica> [hours]}"
HOURS="${3:-11}"
T0=$(date +%s)
D="$RUNS/$SYS/rep$REP"

[ -f "$D/topol.top" ] || [ -f "$(build_dir "$SYS" "$REP")/solv_ions.gro" ] \
  || die "no build for $SYS rep$REP -- run 01_prepare.sh $SYS --replica $REP first"
if [ -f "$D/prod.gro" ]; then log "$SYS rep$REP already finished"; exit 0; fi

if [ ! -f "$D/npt_free.gro" ]; then
  log "== equilibration $SYS rep$REP =="
  "$REPO/scripts/02_equilibrate.sh" "$SYS" "$REP"
else
  log "equilibration already done, skipping"
fi

# Remaining hours for production, keeping a 15 min safety margin
LEFT=$(awk -v h="$HOURS" -v t0="$T0" -v now="$(date +%s)" \
  'BEGIN{r = h - (now - t0) / 3600 - 0.25; if (r < 0.25) r = 0.25; printf "%.2f", r}')
log "== production $SYS rep$REP, -maxh $LEFT =="
"$REPO/scripts/03_production.sh" "$SYS" "$REP" "$LEFT"

pyrun "$REPO/scripts/_registry.py" status --runs "$RUNS" || true
