#!/usr/bin/env bash
# MM/PBSA (M3.5) for one replica of an RNA-bound system.
#
#   PROT_GRP=<n> RNA_GRP=<n> ./run_mmpbsa.sh <system> <replica> [indi]
#
# The checks must be made in a fixed order, starting with the per-residue
# decomposition. If its hot spots do not recover the known operator-contact
# residues of MS2 CP, the setup is wrong and no binding energy computed from it
# is meaningful.

set -euo pipefail
source "$(dirname "$0")/../scripts/lib.sh"
[ "$BRUSAFE_OS" != windows ] || die "gmx_MMPBSA (AmberTools) does not run on native Windows -- use Linux, WSL or macOS for M3.5"
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
CHECKS, IN THIS ORDER:
 1. FINAL_DECOMP.dat -- do the hot spots recover the known operator
    contact residues? If not, stop and correct the setup.
 2. Repeat at indi = 1, 2 and 4 and report the spread as part of the
    error bar. A conclusion that changes between 2 and 4 is not robust.
 3. Repeat on two disjoint frame windows.
 4. Report ddG between systems only; absolute dG values are not quoted.
 5. ddG = -RT ln(ratio); RT = 0.616 kcal/mol at 310 K, so one order
    of magnitude in K_d is about 1.4 kcal/mol.
------------------------------------------------------------------
NOTE
