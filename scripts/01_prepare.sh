#!/usr/bin/env bash
# Build a solvated, ionised, index-tagged system ready for equilibration.
#
#   ./01_prepare.sh <system_name> [--source pdb2gmx|charmm-gui]
#
# Input is expected at  input/<system_name>/  containing either
#   pdb2gmx    : <system_name>.pdb   (cleaned, correct chains, no waters/alt-locs)
#   charmm-gui : the unzipped CHARMM-GUI Solution Builder gromacs/ output
#
# For RNA-containing systems the CHARMM-GUI route is STRONGLY preferred.
# pdb2gmx residue mapping for RNA under the charmm36 port is fragile and the
# Mg2+ placement question below is handled properly by the builder.
set -euo pipefail
source "$(dirname "$0")/lib.sh"
need gmx

SYS="${1:?usage: 01_prepare.sh <system_name> [--source pdb2gmx|charmm-gui]}"
SOURCE="pdb2gmx"
[ "${2:-}" = "--source" ] && SOURCE="${3:?--source needs a value}"

HAS_RNA="$(sysfield "$SYS" has_rna)"
MG="$(sysfield "$SYS" mg_count)"
IN="$REPO/input/$SYS"
OUT="$RUNS/$SYS/build"
mkdir -p "$OUT"; cd "$OUT"
log "preparing $SYS (rna=$HAS_RNA mg=$MG source=$SOURCE)"

if [ "$SOURCE" = "charmm-gui" ]; then
  [ -d "$IN/gromacs" ] || die "expected $IN/gromacs from CHARMM-GUI"
  cp "$IN"/gromacs/step3_input.gro solv_ions.gro
  cp "$IN"/gromacs/topol.top .
  cp -r "$IN"/gromacs/toppar . 2>/dev/null || true
  log "using CHARMM-GUI build as-is (box, salt and Mg2+ set in the builder)"
else
  [ -f "$IN/$SYS.pdb" ] || die "expected $IN/$SYS.pdb"

  # 15 = CHARMM36 in the standard pdb2gmx menu ordering; -ter interactive so you
  # consciously choose termini rather than accepting a default you never saw.
  gmx pdb2gmx -f "$IN/$SYS.pdb" -o proc.gro -p topol.top -i posre.itp \
              -water tip3p -ff charmm36-jul2022 -ignh -ter

  gmx editconf -f proc.gro -o box.gro -c -d 1.2 -bt dodecahedron
  gmx solvate  -cp box.gro -cs spc216.gro -o solv.gro -p topol.top

  gmx grompp -f "$MDP/em.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 1

  if [ "$HAS_RNA" = "yes" ] && [ "$MG" -gt 0 ]; then
    # ---------------------------------------------------------------------
    # Mg2+ WARNING -- read before trusting anything downstream.
    # Mg2+ water exchange is a microsecond process. On a 250 ns trajectory the
    # ions never equilibrate: they stay essentially wherever genion drops them.
    # Random placement therefore adds noise, not realism.
    # Preferred order:
    #   1. crystallographic Mg2+ sites from the cocrystal, placed explicitly
    #   2. CHARMM-GUI placement with its ion-distribution options
    #   3. genion (this branch) -- acceptable only if you document it and run
    #      04_postprocess.sh's mg_check, which flags Mg2+ that wandered into
    #      the binding interface.
    # ---------------------------------------------------------------------
    log "WARNING: genion Mg2+ placement is arbitrary and will not equilibrate"
    printf 'SOL\n' | gmx genion -s ions.tpr -o mg.gro -p topol.top \
                                -pname MG -pq 2 -np "$MG"
    gmx grompp -f "$MDP/em.mdp" -c mg.gro -p topol.top -o ions2.tpr -maxwarn 1
    printf 'SOL\n' | gmx genion -s ions2.tpr -o solv_ions.gro -p topol.top \
                                -pname POT -nname CLA -conc 0.15 -neutral
  else
    printf 'SOL\n' | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top \
                                -pname POT -nname CLA -conc 0.15 -neutral
  fi
fi

# ---- index groups -----------------------------------------------------------
# TWO sources, concatenated. Both are needed and neither is optional:
#
#   default.ndx  from make_ndx -- System, Protein, Backbone, C-alpha, RNA, ...
#                04_postprocess.sh asks for these BY NAME (Backbone, C-alpha,
#                RNA). `gmx select -on` does not emit them.
#   custom.ndx   from gmx select -- SOLU / SOLV. Every mdp couples to these,
#                so one mdp set works for apo, RNA-bound and RNA-only systems.
#
# .ndx is plain text ([ name ] followed by atom numbers), so `cat` is a valid
# merge. Group names must stay unique across the two files.
EXCL="$SOLVENT_RESNAMES"
printf 'q\n' | gmx make_ndx -f solv_ions.gro -o default.ndx >/dev/null 2>&1
gmx select -s solv_ions.gro -on custom.ndx \
  -select "\"SOLU\" not resname $EXCL; \"SOLV\" resname $EXCL"
cat default.ndx custom.ndx > index.ndx
rm -f default.ndx custom.ndx

log "index groups: $(grep -o '\[ [^]]* \]' index.ndx | tr -d '[]' | tr -s ' ' | paste -sd' ' -)"

# Fail loudly here rather than three days later inside 04_postprocess.sh.
for g in Backbone C-alpha SOLU SOLV; do
  grep -q "\[ $g \]" index.ndx || die "index.ndx has no '$g' group -- 04_postprocess.sh will fail"
done
if [ "$HAS_RNA" = "yes" ]; then
  grep -q '\[ RNA \]' index.ndx || log "WARNING: no 'RNA' group. Add one before postprocessing:
      gmx make_ndx -f $OUT/solv_ions.gro -n $OUT/index.ndx -o $OUT/index.ndx"
fi

log "built $OUT/solv_ions.gro and $OUT/index.ndx"
log "next: 02_equilibrate.sh $SYS <replica>"
