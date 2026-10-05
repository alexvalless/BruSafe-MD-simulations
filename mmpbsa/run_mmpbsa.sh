#!/usr/bin/env bash

# Gate order matters: decomposition first. If the per-residue hot spots do not
# recover the known operator-contact residues of MS2 CP, the setup is wrong
# and no binding energy from it means anything. Do not skip ahead.

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

# Frame window from the registry: drop skip_ns, then thin to ~MMPBSA_FRAMES
# frames per replica. Error bars come from the replica spread, so more frames
# from one replica mostly buy CPU time, not information.
PMDP="$D/prod.mdp"
FRAME_PS="$(awk -v n="$(mdpval "$PMDP" nstxout-compressed)" -v dt="$(mdpval "$PMDP" dt)" 'BEGIN{print n * dt}')"
read -r START INTERVAL < <(awk -v skip="$(sysfield "$SYS" skip_ns)" -v ns="$(sysfield "$SYS" ns)" \
  -v fps="$FRAME_PS" -v want="${MMPBSA_FRAMES:-250}" 'BEGIN{
    s = int(skip * 1000 / fps + 0.5) + 1; n = (ns - skip) * 1000 / fps
    i = int(n / want); if (i < 1) i = 1; print s, i }')
log "frames: start=$START interval=$INTERVAL (frame = $FRAME_PS ps)"

sed "s/^  indi .*/  indi                  = $INDI/" "$REPO/mmpbsa/mmpbsa_M3.5.in" > mmpbsa.in
sed -i -e "s/^  startframe .*/  startframe            = $START/" \
       -e "s/^  endframe .*/  endframe              = 9999999/" \
       -e "s/^  interval .*/  interval              = $INTERVAL/" mmpbsa.in
sed -i "s/^  sys_name .*/  sys_name              = \"${SYS}_rep${REP}\"/" mmpbsa.in

# Group numbers looked up by name unless given explicitly
PROT_GRP="${PROT_GRP:-$(ndxgroup "$D/index.ndx" Protein)}"
RNA_GRP="${RNA_GRP:-$(ndxgroup "$D/index.ndx" RNA)}"
[ -n "$PROT_GRP" ] && [ -n "$RNA_GRP" ] || die "no Protein / RNA group in $D/index.ndx"
log "groups: Protein=$PROT_GRP RNA=$RNA_GRP"

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
