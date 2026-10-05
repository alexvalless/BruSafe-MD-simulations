#!/usr/bin/env bash
# Record the exact software/hardware state of THIS machine. Run once per box,
# commit the output. If two machines disagree, stop and fix it before running
# anything -- mixed GROMACS versions across replicas invalidates the campaign.
set -euo pipefail
source "$(dirname "$0")/lib.sh"

OUT="$REPO/docs/env_$(hostname).txt"
{
  echo "host        : $(hostname)"
  echo "date        : $(now_iso)"
  echo "os          : $BRUSAFE_OS ($(uname -s) $(uname -r))"
  echo "cpu         : $( (grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | xargs) || sysctl -n machdep.cpu.brand_string 2>/dev/null || echo unknown)"
  echo "cores       : $(ncpu)"
  echo "python      : ${PY:-NOT FOUND} $( [ -n "$PY" ] && pyrun -c 'import sys; print(sys.version.split()[0])' )"
  echo "numpy       : $( [ -n "$PY" ] && pyrun -c 'import numpy; print(numpy.__version__)' 2>/dev/null || echo NOT FOUND)"
  echo
  echo "--- gromacs ---"
  echo "binary: $GMX_BIN"
  gmx --version 2>/dev/null | tr -d '\r' | sed -n '1,25p' || echo "gmx NOT FOUND"
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
