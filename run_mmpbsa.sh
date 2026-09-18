#!/usr/bin/env bash
# M3.5 -- MM/PBSA on one RNA-bound system.
#   ./run_mmpbsa.sh <system_name> <replica> [indi]
#
# Gate order matters: decomposition FIRST. If the per-residue hot spots do not
# recover the known operator-contact residues of MS2 CP, the setup is wrong
# and no binding energy from it means anything. Do not skip ahead.
set -euo pipefail
source "$(dirname "$0")/../scripts/lib.sh"
need gmx_MMPBSA

SYS="${1:?usage: run_mmpbsa.sh <system> <replica> [indi]}"
REP="${2:?usage: run_mmpbsa.sh <system> <replica> [indi]}"
INDI="${3:-4.0}"
[ "$(sysfield "$SYS" has_rna)" = "yes" ] || die "$SYS has no RNA -- nothing to bind"

D="$RUNS/$SYS/rep$REP"
OUT="$D/mmpbsa_indi${INDI}"
mkdir -p "$OUT"; cd "$OUT"
[ -f "$D/analysis/clean.xtc" ] || die "run 04_postprocess.sh $SYS $REP first"

log "stripping solvent (implicit-solvent method -- water must not be present)"
printf 'SOLU\n' | gmx trjconv -s "$D/prod.tpr" -f "$D/analysis/clean.xtc" \
                              -n "$D/index.ndx" -o complex_dry.xtc

sed "s/^  indi .*/  indi                  = $INDI/" "$REPO/mmpbsa/mmpbsa_M3.5.in" > mmpbsa.in
sed -i "s/^  sys_name .*/  sys_name              = \"${SYS}_rep${REP}\"/" mmpbsa.in

log "running gmx_MMPBSA (indi=$INDI) -- PB primary, GB cross-check"
mpirun -np "${MMPBSA_NP:-8}" gmx_MMPBSA MPI -O -i mmpbsa.in \
  -cs "$D/prod.tpr" -ci "$D/index.ndx" -cg "${PROT_GRP:?set PROT_GRP}" "${RNA_GRP:?set RNA_GRP}" \
  -ct complex_dry.xtc -cp "$D/topol.top" \
  -o FINAL_RESULTS.dat -eo FINAL_RESULTS.csv \
  -do FINAL_DECOMP.dat -deo FINAL_DECOMP.csv

cat <<'NOTE'

------------------------------------------------------------------
CHECK IN THIS ORDER:
 1. FINAL_DECOMP.dat -- do the hot spots recover the known operator
    contact residues? If not, STOP. Fix the setup.
 2. Repeat at indi = 1, 2 and 4. Report the spread as part of the
    error bar. A conclusion that flips between 2 and 4 is not a
    conclusion.
 3. Repeat on two disjoint frame windows.
 4. Report ddG between systems only. Never quote an absolute dG.
 5. ddG = -RT ln(ratio); RT = 0.616 kcal/mol at 310 K, so one order
    of magnitude in K_d is about 1.4 kcal/mol.
------------------------------------------------------------------
NOTE
