#!/usr/bin/env bash
# EM -> NVT(posres) -> NPT(posres) -> NPT(free), for ONE replica.
#
#   ./02_equilibrate.sh <system_name> <replica_number>
#
# The velocity seed is generated per replica here. This is what makes replicas
# statistically independent. Branching three replicas off one equilibrated
# frame gives you one run wearing a trenchcoat -- do not do it.
set -euo pipefail
source "$(dirname "$0")/lib.sh"
need gmx

SYS="${1:?usage: 02_equilibrate.sh <system> <replica>}"
REP="${2:?usage: 02_equilibrate.sh <system> <replica>}"
BUILD="$RUNS/$SYS/build"
D="$RUNS/$SYS/rep$REP"
[ -f "$BUILD/solv_ions.gro" ] || die "run 01_prepare.sh $SYS first"
mkdir -p "$D"; cd "$D"
guard_existing "$D"

SEED="$(replica_seed "$SYS" "$REP")"
log "$SYS rep$REP  seed=$SEED"
echo "$SEED" > seed.txt

cp "$BUILD/topol.top" "$BUILD/index.ndx" .
cp -r "$BUILD/toppar" . 2>/dev/null || true
cp "$BUILD"/*.itp . 2>/dev/null || true
sed "s/REPLACE_SEED/$SEED/" "$MDP/nvt.mdp" > nvt.mdp

log "energy minimisation"
gmx grompp -f "$MDP/em.mdp" -c "$BUILD/solv_ions.gro" -r "$BUILD/solv_ions.gro" \
           -p topol.top -n index.ndx -o em.tpr -maxwarn 1
gmx mdrun -deffnm em -ntmpi 1

log "NVT 200 ps (restrained)"
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -deffnm nvt $MDRUN_FLAGS

log "NPT 1 ns (restrained)"
gmx grompp -f "$MDP/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt \
           -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -deffnm npt $MDRUN_FLAGS

log "NPT 5 ns (unrestrained)"
gmx grompp -f "$MDP/npt_free.mdp" -c npt.gro -t npt.cpt \
           -p topol.top -n index.ndx -o npt_free.tpr
gmx mdrun -deffnm npt_free $MDRUN_FLAGS

# ---- equilibration gate --------------------------------------------------
log "gate: density / temperature / pressure"
printf 'Density\n\n'     | gmx energy -f npt_free.edr -o density.xvg     >/dev/null 2>&1
printf 'Temperature\n\n' | gmx energy -f npt_free.edr -o temperature.xvg >/dev/null 2>&1
printf 'Pressure\n\n'    | gmx energy -f npt_free.edr -o pressure.xvg    >/dev/null 2>&1
python3 "$REPO/analysis/convergence.py" gate \
  --density density.xvg --temperature temperature.xvg --ref-t 310 \
  || die "equilibration gate FAILED for $SYS rep$REP -- do not start production"

log "equilibration complete: $D/npt_free.gro"
