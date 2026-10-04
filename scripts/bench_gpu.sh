#!/usr/bin/env bash
# Throughput sweep on one GPU. Run ONCE, on S1, before the campaign starts.
# Output decides MDRUN_FLAGS for every subsequent run.
# usage: bench_gpu.sh <path/to/prod.tpr>
set -euo pipefail
source "$(dirname "$0")/lib.sh"
need gmx

TPR="${1:?usage: bench_gpu.sh <prod.tpr>}"
TPR="$(readlink -f "$TPR")"
NSTEPS=60000
BDIR="$REPO/bench_$(hostname)"
mkdir -p "$BDIR"; cd "$BDIR"
OUT="results.txt"; : > "$OUT"

nvidia-smi --query-gpu=name,memory.total,power.limit --format=csv,noheader | tee -a "$OUT"
echo "cores: $(nproc)" | tee -a "$OUT"; echo | tee -a "$OUT"

run() {
  local label="$1" omp="$2"; shift 2
  local d="cfg_${label}_omp${omp}"
  rm -rf "$d"; mkdir -p "$d"
  if ( cd "$d" && gmx mdrun -s "$TPR" -deffnm b -ntmpi 1 -ntomp "$omp" -pin on \
        -nsteps "$NSTEPS" -resethway -noconfout -nstlist 200 "$@" >/dev/null 2>&1 ); then
    printf '%-16s omp=%-3s %9s ns/day\n' "$label" "$omp" \
      "$(awk '/Performance:/{print $2}' "$d/b.log")" | tee -a "$OUT"
  else
    printf '%-16s omp=%-3s %9s\n' "$label" "$omp" FAILED | tee -a "$OUT"
  fi
}

echo "--- offload configurations ---" | tee -a "$OUT"
for omp in 4 6 8 12; do
  [ "$omp" -le "$(nproc)" ] || continue
  run nb_pme        "$omp" -nb gpu -pme gpu
  run nb_pme_bonded "$omp" -nb gpu -pme gpu -bonded gpu
  run gpu_resident  "$omp" -nb gpu -pme gpu -bonded gpu -update gpu
done

echo | tee -a "$OUT"
echo "--- 2 concurrent runs, same GPU (aggregate is what matters) ---" | tee -a "$OUT"
H=$(( $(nproc) / 2 )); [ "$H" -lt 2 ] && H=2
for i in 0 1; do
  d="conc_$i"; rm -rf "$d"; mkdir -p "$d"
  ( cd "$d" && gmx mdrun -s "$TPR" -deffnm b -ntmpi 1 -ntomp "$H" -pin on \
      -pinoffset $(( i * H )) -nb gpu -pme gpu -bonded gpu -update gpu \
      -nsteps "$NSTEPS" -resethway -noconfout -nstlist 200 >/dev/null 2>&1 ) &
done
wait
T=0
for i in 0 1; do
  p=$(awk '/Performance:/{print $2}' "conc_$i/b.log")
  printf '  run %d: %9s ns/day\n' "$i" "$p" | tee -a "$OUT"
  T=$(python3 -c "print($T + $p)")
done
printf '  AGGREGATE: %9.1f ns/day\n' "$T" | tee -a "$OUT"

cat >> "$OUT" <<'NOTE'

NEXT STEP: export the winning configuration and never change it.
  export BRUSAFE_MDRUN_FLAGS="-nb gpu -pme gpu -bonded gpu -update gpu -ntmpi 1 -ntomp <N> -pin on"
Record it in docs/DECISIONS.md with the date.
NOTE
cat "$OUT"
