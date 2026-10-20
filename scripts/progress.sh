#!/usr/bin/env bash
# One-screen progress report for THIS machine. Safe to run while mdrun is going.
#
#   ./scripts/progress.sh          # print once
#   ./scripts/progress.sh watch    # refresh every 30 s (Ctrl+C to leave; the run is not affected)
#
# Shows: every replica under runs/ that has started, for each running
# production how big prod.xtc is and how long ago it last grew (a stuck run
# stops growing), the simulated time at the last checkpoint (written every
# 15 min), the GPU load, and the tail of the slot log if there is one.
# On Windows mdrun keeps prod.log locked, so the checkpoint is the only
# place to read the exact simulated time.
set -uo pipefail
source "$(dirname "$0")/lib.sh"

mtime() { stat -c %Y "$1" 2>/dev/null || stat -f %m "$1" 2>/dev/null || echo 0; }

report() {
  local now; now="$(date +%s)"
  echo "=== $(hostname)   $(date '+%Y-%m-%d %H:%M:%S') ==="
  echo
  pyrun "$REPO/scripts/_registry.py" status --runs "$RUNS" 2>&1 \
    | grep -vE 'NOT STARTED|never launched|^  [A-Za-z0-9_]+/rep' || true
  echo

  local d any=0 sys rep ns target cpt t xtc age mb
  for d in "$RUNS"/*/rep*; do
    [ -d "$d" ] || continue
    [ -f "$d/prod.gro" ] && continue
    [ -f "$d/prod.log" ] || continue
    any=1
    sys="$(basename "$(dirname "$d")")"; rep="${d##*rep}"
    target="$(sysfield "$sys" ns 2>/dev/null || echo '?')"
    xtc="$d/prod.xtc"
    if [ -f "$xtc" ]; then
      age=$(( now - $(mtime "$xtc") ))
      mb=$(( $(wc -c < "$xtc" | tr -d ' ') / 1048576 ))
    else
      age="-"; mb=0
    fi
    ns="-"
    if [ -f "$d/prod.cpt" ]; then
      t="$( (gmx dump -cp "$d/prod.cpt" 2>/dev/null || true) | awk '/^t = /{v=$3} END{if (v != "") print v}' | tr -d '\r')"
      [ -n "$t" ] && ns="$(awk -v t="$t" 'BEGIN{printf "%.2f", t / 1000}')"
    fi
    printf '  %-18s rep%-2s  checkpoint %6s / %s ns   prod.xtc %5s MB, last grew %ss ago\n' \
           "$sys" "$rep" "$ns" "$target" "$mb" "$age"
    if [ "$age" != "-" ] && [ "$age" -gt 600 ]; then
      echo "      ^ the trajectory has not grown for $((age / 60)) min: the run may be stuck or stopped"
    fi
  done
  [ "$any" = 1 ] || echo "  (no production running here right now)"
  echo

  if command -v nvidia-smi >/dev/null 2>&1; then
    echo "GPU: $(nvidia-smi --query-gpu=utilization.gpu,temperature.gpu,power.draw --format=csv,noheader | tr -d '\r')"
  fi
  local f
  for f in "$REPO"/slot_*.log; do
    [ -f "$f" ] || continue
    echo
    echo "$(basename "$f"):"
    grep -E '^\[' "$f" 2>/dev/null | tail -4 | cut -c1-140 | sed 's/^/  /' || true
  done
}

if [ "${1:-}" = watch ]; then
  while true; do
    printf '\033c'
    report
    echo; echo "(refreshing every 30 s, Ctrl+C to leave)"
    sleep 30
  done
else
  report
fi
