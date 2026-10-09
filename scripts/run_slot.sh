#!/usr/bin/env bash
# Run every job assigned to ONE machine slot in docs/PLAN.txt, one after another,
# unattended: build the system if needed, equilibrate, produce, post-process.
#
#   ./scripts/run_slot.sh <slot_number> [hours_budget]
#
# Example: this PC is "slot 3" and the night is 16 h long:
#   ./scripts/run_slot.sh 3 16
#
# docs/PLAN.txt comes from:
#   python3 scripts/_registry.py plan --machines 15 --tier 2 --nsday 580 > docs/PLAN.txt
# and lists the jobs as "gpu<slot><TAB>system<TAB>replica<TAB>ns<TAB>consumer".
#
# Unattended builds only work for systems with a CHARMM-GUI build in
# input/<system>/gromacs/ (pdb2gmx asks about termini interactively).
# Safe to rerun after a crash or reboot: finished work is skipped and an
# interrupted production resumes from its checkpoint.
set -euo pipefail
source "$(dirname "$0")/lib.sh"
need gmx

SLOT="${1:?usage: run_slot.sh <slot_number> [hours_budget]}"
HOURS="${2:-14}"
PLAN="$REPO/docs/PLAN.txt"
[ -f "$PLAN" ] || die "no $PLAN -- generate it with scripts/_registry.py plan"

T0=$(date +%s)
JOBS="$(awk -F'\t' -v s="gpu$SLOT" '$1 == s {print $2 "\t" $3}' "$PLAN" | tr -d '\r')"
[ -n "$JOBS" ] || die "slot $SLOT has no jobs in docs/PLAN.txt"
log "slot $SLOT on $(hostname): $(printf '%s\n' "$JOBS" | tr '\t' ' ' | paste -sd',' -)"

left_hours() {
  awk -v h="$HOURS" -v t0="$T0" -v now="$(date +%s)" \
    'BEGIN{printf "%.2f", h - (now - t0) / 3600}'
}

FAILED=""
while IFS=$'\t' read -r SYS REP; do
  [ -n "$SYS" ] || continue
  D="$RUNS/$SYS/rep$REP"
  LEFT="$(left_hours)"
  if awk -v l="$LEFT" 'BEGIN{exit !(l < 0.5)}'; then
    log "less than 30 min of budget left, not starting $SYS rep$REP"; FAILED="$FAILED $SYS/rep$REP"; continue
  fi

  if [ ! -f "$D/prod.gro" ]; then
    # a replica already started keeps the build it started with; a new one
    # gets its own ion placement (build_rep<N>), see lib.sh build_dir
    if [ ! -f "$D/topol.top" ] && [ ! -f "$RUNS/$SYS/build_rep$REP/solv_ions.gro" ]; then
      [ -d "$REPO/input/$SYS/gromacs" ] || { log "NO CHARMM-GUI build for $SYS (input/$SYS/gromacs) -- skipped"; FAILED="$FAILED $SYS/rep$REP"; continue; }
      log "== building $SYS for replica $REP =="
      mkdir -p "$RUNS"
      "$REPO/scripts/01_prepare.sh" "$SYS" --source charmm-gui --replica "$REP" > "$RUNS/prepare_${SYS}_rep${REP}.log" 2>&1 \
        || { log "prepare FAILED for $SYS rep$REP (see $RUNS/prepare_${SYS}_rep${REP}.log)"; FAILED="$FAILED $SYS/rep$REP"; continue; }
    fi
    "$REPO/scripts/run_night.sh" "$SYS" "$REP" "$LEFT" \
      || { log "run FAILED for $SYS rep$REP"; FAILED="$FAILED $SYS/rep$REP"; continue; }
  fi

  if [ ! -f "$D/prod.gro" ]; then
    log "$SYS rep$REP did not finish within the budget; rerun this command to resume"
    FAILED="$FAILED $SYS/rep$REP"; break
  fi
  if [ ! -f "$D/analysis/clean.xtc" ]; then
    log "== post-processing $SYS rep$REP =="
    "$REPO/scripts/04_postprocess.sh" "$SYS" "$REP" > "$D/postprocess.log" 2>&1 \
      || log "postprocess had errors (see $D/postprocess.log); the trajectory itself is intact"
  fi
done <<< "$JOBS"

# pack whatever finished (and copy to OneDrive if it is set up); never fatal
"$REPO/scripts/pack_results.sh" > "$RUNS/pack_slot${SLOT}.log" 2>&1 \
  && log "results packed (see pack/ and $RUNS/pack_slot${SLOT}.log)" \
  || log "packing failed (see $RUNS/pack_slot${SLOT}.log); run ./scripts/pack_results.sh by hand"

log "slot $SLOT done in $(awk -v t0="$T0" -v now="$(date +%s)" 'BEGIN{printf "%.1f", (now - t0) / 3600}') h"
if [ -n "$FAILED" ]; then log "NOT completed:$FAILED"; exit 1; fi
log "all jobs of slot $SLOT completed"
