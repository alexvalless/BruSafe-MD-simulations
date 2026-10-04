#!/usr/bin/env bash
# Production MD for ONE replica. Restartable: rerun the same command after a
# crash or reboot and it picks up from the last checkpoint.
#
#   ./03_production.sh <system_name> <replica_number> [max_hours]
set -euo pipefail
source "$(dirname "$0")/lib.sh"
need gmx

SYS="${1:?usage: 03_production.sh <system> <replica> [max_hours]}"
REP="${2:?usage: 03_production.sh <system> <replica> [max_hours]}"
MAXH="${3:-24}"
D="$RUNS/$SYS/rep$REP"
cd "$D" 2>/dev/null || die "no run dir $D -- run 02_equilibrate.sh first"
[ -f npt_free.gro ] || die "equilibration incomplete in $D"

if [ -f prod.gro ]; then log "$SYS rep$REP already finished"; exit 0; fi

if [ ! -f prod.tpr ]; then
  NS="$(sysfield "$SYS" ns)"
  STEPS=$(( NS * 500000 ))            # ns -> steps at dt = 2 fs
  log "grompp for $NS ns ($STEPS steps)"
  sed "s/^nsteps .*/nsteps                  = $STEPS/" "$MDP/prod.mdp" > prod.mdp
  gmx grompp -f prod.mdp -c npt_free.gro -t npt_free.cpt \
             -p topol.top -n index.ndx -o prod.tpr
  # provenance: what actually ran, alongside the trajectory
  { echo "host      : $(hostname)"; echo "started   : $(date -Is)";
    echo "gmx       : $(gmx --version 2>/dev/null | head -1)";
    echo "flags     : $MDRUN_FLAGS"; echo "seed      : $(cat seed.txt)";
    echo "target_ns : $NS"; } > PROVENANCE.txt
fi

if [ -f prod.cpt ]; then
  log "resuming $SYS rep$REP from checkpoint"
  gmx mdrun -deffnm prod -cpi prod.cpt -append $MDRUN_FLAGS -maxh "$MAXH" -cpt 15
else
  log "starting $SYS rep$REP"
  gmx mdrun -deffnm prod $MDRUN_FLAGS -maxh "$MAXH" -cpt 15
fi

if [ -f prod.gro ]; then
  echo "finished  : $(date -Is)" >> PROVENANCE.txt
  log "COMPLETE $SYS rep$REP"
else
  log "hit wall clock -- rerun the same command to continue"
fi
