#!/usr/bin/env bash
# Record the exact software/hardware state of THIS machine. Run once per box,
# commit the output. If two machines disagree, stop and fix it before running
# anything -- mixed GROMACS versions across replicas invalidates the campaign.
set -euo pipefail
source "$(dirname "$0")/lib.sh"

OUT="$REPO/docs/env_$(hostname).txt"
{
  echo "host        : $(hostname)"
  echo "date        : $(date -Is)"
  echo "kernel      : $(uname -r)"
  echo "cpu         : $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)"
  echo "cores       : $(nproc)"
  echo
  echo "--- gromacs ---"
  gmx --version 2>/dev/null | sed -n '1,25p' || echo "gmx NOT FOUND"
  echo
  echo "--- gpu ---"
  nvidia-smi --query-gpu=name,driver_version,memory.total,power.limit \
             --format=csv 2>/dev/null || echo "nvidia-smi NOT FOUND"
  echo
  echo "--- force field ---"
  echo "GMXLIB=${GMXLIB:-<unset>}"
  ls -d "${GMXLIB:-.}"/charmm36*.ff 2>/dev/null || echo "charmm36 ff NOT in GMXLIB"
  echo
  echo "--- locked mdrun flags ---"
  echo "$MDRUN_FLAGS"
} | tee "$OUT"

log "wrote $OUT -- commit this file"
