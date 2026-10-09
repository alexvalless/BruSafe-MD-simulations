#!/usr/bin/env bash
# Build a solvated, ionised, index-tagged system ready for equilibration.
#
#   ./01_prepare.sh <system_name> [--source pdb2gmx|charmm-gui]
#
# Input is expected at  input/<system_name>/  containing either
#   pdb2gmx    : <system_name>.pdb   (cleaned, correct chains, no waters/alt-locs)
#   charmm-gui : the unzipped CHARMM-GUI Solution Builder gromacs/ output
#
# For RNA-containing systems the CHARMM-GUI route is STRONGLY preferred:
# pdb2gmx residue mapping for RNA under the charmm36 port is fragile. Only the
# CHARMM-GUI topology is used; box, water and ions are rebuilt here the same
# way for both routes.
set -euo pipefail
source "$(dirname "$0")/lib.sh"
need gmx

SYS="${1:?usage: 01_prepare.sh <system_name> [--source pdb2gmx|charmm-gui]}"
SOURCE="pdb2gmx"
[ "${2:-}" = "--source" ] && SOURCE="${3:?--source needs a value}"

# Box padding (nm) from the registry's box_nm column: 1.0 keeps periodic images
# 2.0 nm apart, above the 1.2 nm cutoff; the free hairpin (S8) gets 1.2.
BOX_PAD="${BRUSAFE_BOX_PAD:-$(sysfield "$SYS" box_nm)}"
FF="${BRUSAFE_FF:-charmm36-jul2022}"
HAS_RNA="$(sysfield "$SYS" has_rna)"
MG="$(sysfield "$SYS" mg_count)"
IN="$REPO/input/$SYS"
OUT="$RUNS/$SYS/build"
mkdir -p "$OUT"; cd "$OUT"
log "preparing $SYS (rna=$HAS_RNA mg=$MG source=$SOURCE pad=${BOX_PAD}nm)"

if [ "$SOURCE" = "charmm-gui" ]; then
  # ---------------------------------------------------------------------
  # CHARMM-GUI supplies the TOPOLOGY (CHARMM36m parameters, termini, rebuilt
  # atoms). The box is rebuilt here with GROMACS: CHARMM-GUI's octahedral
  # GROMACS output was found to place hundreds of atoms on top of their own
  # periodic images (energy minimisation diverges), and rebuilding also gives
  # every system the same rhombic-dodecahedron protocol as the pdb2gmx route.
  # BRUSAFE_KEEP_CGUI_BOX=1 keeps CHARMM-GUI's box (only for a RECTANGULAR one).
  # ---------------------------------------------------------------------
  G="$IN/gromacs"
  [ -f "$G/step3_input.gro" ] && [ -f "$G/topol.top" ] || die "expected $G from CHARMM-GUI"
  rm -rf toppar; cp -r "$G/toppar" .
  if [ "${BRUSAFE_KEEP_CGUI_BOX:-0}" = "1" ]; then
    cp "$G/step3_input.gro" solv_ions.gro
    cp "$G/topol.top" topol.top
    log "using CHARMM-GUI box as-is (BRUSAFE_KEEP_CGUI_BOX=1)"
  else
    log "CHARMM-GUI topology; rebuilding box, water and ions with GROMACS"
    gmx select -s "$G/step3_input.gro" -on solute.ndx \
               -select "\"solute\" not resname $SOLVENT_RESNAMES"
    printf '0\n' | gmx editconf -f "$G/step3_input.gro" -n solute.ndx -o proc.gro
    # same topology, minus every solvent / ion molecule line
    awk -v drop=" $SOLVENT_RESNAMES " '
      /^[ \t]*\[/ { insec = ($0 ~ /molecules/) }
      insec && NF == 2 && $1 !~ /^;/ && index(drop, " " $1 " ") { next }
      { print }' "$G/topol.top" > topol.top
  fi
else
  [ -f "$IN/$SYS.pdb" ] || die "expected $IN/$SYS.pdb"

  check_ff "$FF"

  # -ter asks questions on the terminal. Git Bash's default window (mintty)
  # does not pass keyboard input to native Windows programs; winpty fixes it.
  WINPTY=""
  if [ "$BRUSAFE_OS" = windows ] && [ -t 0 ] && command -v winpty >/dev/null 2>&1; then
    WINPTY="winpty"
  fi

  # 15 = CHARMM36 in the standard pdb2gmx menu ordering; -ter interactive so you
  # consciously choose termini rather than accepting a default you never saw.
  $WINPTY "$GMX_BIN" pdb2gmx -f "$IN/$SYS.pdb" -o proc.gro -p topol.top -i posre.itp \
              -water tip3p -ff "$FF" -ignh -ter
fi

if [ ! -f solv_ions.gro ]; then
  gmx editconf -f proc.gro -o box.gro -c -d "$BOX_PAD" -bt dodecahedron
  gmx solvate  -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
  WATER=SOL
  if [ "$SOURCE" = "charmm-gui" ]; then
    # CHARMM-GUI names its water TIP3 / OH2 H1 H2; rename the solvate output
    awk 'NR > 2 && substr($0, 6, 5) == "SOL  " {
           a = substr($0, 11, 5); gsub(/ /, "", a)
           n = (a == "OW") ? "OH2" : (a == "HW1") ? "H1" : (a == "HW2") ? "H2" : a
           $0 = substr($0, 1, 5) "TIP3 " sprintf("%5s", n) substr($0, 16) }
         { print }' solv.gro > solv.tmp && mv solv.tmp solv.gro
    sed -i.bak -E 's/^SOL([[:space:]]+[0-9]+[[:space:]]*)$/TIP3\1/' topol.top && rm -f topol.top.bak
    WATER=TIP3
  fi

  gmx grompp -f "$MDP/em.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 1

  if [ "$HAS_RNA" = "yes" ] && [ "$MG" -gt 0 ]; then
    # ---------------------------------------------------------------------
    # Mg2+ WARNING -- read before trusting anything downstream.
    # Mg2+ water exchange is a microsecond process. On a 30-100 ns trajectory the
    # ions never equilibrate: they stay essentially wherever genion drops them.
    # Random placement therefore adds noise, not realism. Crystallographic
    # sites would be better; 1ZDH has none, so this is documented in
    # docs/DECISIONS.md and 04_postprocess.sh runs the interface check
    # (mg_interface.txt), which flags Mg2+ sitting in the binding interface.
    # ---------------------------------------------------------------------
    log "WARNING: genion Mg2+ placement is arbitrary and will not equilibrate"
    printf '%s\n' "$WATER" | gmx genion -seed 2026 -s ions.tpr -o mg.gro -p topol.top \
                                -pname MG -pq 2 -np "$MG"
    gmx grompp -f "$MDP/em.mdp" -c mg.gro -p topol.top -o ions2.tpr -maxwarn 1
    printf '%s\n' "$WATER" | gmx genion -seed 2026 -s ions2.tpr -o solv_ions.gro -p topol.top \
                                -pname POT -nname CLA -conc 0.15 -neutral
  else
    printf '%s\n' "$WATER" | gmx genion -seed 2026 -s ions.tpr -o solv_ions.gro -p topol.top \
                                -pname POT -nname CLA -conc 0.15 -neutral
  fi
fi

# ---- hydrogen mass repartitioning --------------------------------------------
# Written into the local topology files (never the GMXLIB force field), so it
# works on any GROMACS version and dt = 4 fs is safe in every mdp. Skips itself
# if the topology is already repartitioned (e.g. CHARMM-GUI's HMR option).
pyrun "$REPO/scripts/hmr_top.py" topol.top

# ---- index groups -----------------------------------------------------------
# TWO sources, concatenated. Both are needed and neither is optional:
#
#   default.ndx  from make_ndx -- System, Protein, Backbone, C-alpha, RNA, ...
#                04_postprocess.sh asks for these BY NAME (Backbone, C-alpha,
#                RNA). `gmx select -on` does not emit them.
#   custom.ndx   from gmx select -- SOLU / SOLV. Every mdp couples to these,
#                so one mdp set works for apo, RNA-bound and RNA-only systems.
#                RNA systems also get RNA, RNA_BB (sugar-phosphate backbone)
#                and RNA_P (phosphorus); 04_postprocess.sh uses RNA_BB / RNA_P
#                in place of Backbone / C-alpha when there is no protein.
#
# .ndx is plain text ([ name ] followed by atom numbers), so `cat` is a valid
# merge. Group names must stay unique across the two files.
EXCL="$SOLVENT_RESNAMES"
printf 'q\n' | gmx make_ndx -f solv_ions.gro -o default.ndx >/dev/null 2>&1
SEL="\"SOLU\" not resname $EXCL; \"SOLV\" resname $EXCL"
if [ "$HAS_RNA" = "yes" ]; then
  R="resname $RNA_RESNAMES"
  grep -q '\[ RNA \]' default.ndx || SEL="$SEL; \"RNA\" $R"
  SEL="$SEL; \"RNA_BB\" $R and name P \"O5'\" \"C5'\" \"C4'\" \"C3'\" \"O3'\""
  SEL="$SEL; \"RNA_P\" $R and name P"
fi
gmx select -s solv_ions.gro -on custom.ndx -select "$SEL"
cat default.ndx custom.ndx > index.ndx
rm -f default.ndx custom.ndx

log "index groups: $(grep -o '\[ [^]]* \]' index.ndx | tr -d '[]' | tr -s ' ' | paste -sd' ' -)"

# Fail loudly here rather than three days later inside 04_postprocess.sh.
NEED="SOLU SOLV"
grep -q '\[ Protein \]' index.ndx && NEED="$NEED Backbone C-alpha"
[ "$HAS_RNA" = "yes" ] && NEED="$NEED RNA RNA_BB RNA_P"
for g in $NEED; do
  grep -q "\[ $g \]" index.ndx || die "index.ndx has no '$g' group -- 04_postprocess.sh will fail"
done

log "built $OUT/solv_ions.gro and $OUT/index.ndx"
log "next: 02_equilibrate.sh $SYS <replica>"
