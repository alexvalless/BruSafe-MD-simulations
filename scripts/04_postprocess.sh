#!/usr/bin/env bash
# PBC cleanup + the standard analysis battery for one replica.
# Skipping the PBC step is the single most common way to produce confidently
# wrong RMSD/MM-PBSA numbers, so it is not optional and not separable.
#
#   ./04_postprocess.sh <system_name> <replica_number>
set -euo pipefail
source "$(dirname "$0")/lib.sh"
need gmx

SYS="${1:?usage: 04_postprocess.sh <system> <replica>}"
REP="${2:?usage: 04_postprocess.sh <system> <replica>}"
HAS_RNA="$(sysfield "$SYS" has_rna)"
SKIP_PS=$(( $(sysfield "$SYS" skip_ns) * 1000 ))
D="$RUNS/$SYS/rep$REP"
A="$D/analysis"
cd "$D"; mkdir -p "$A"
[ -f prod.xtc ] || die "no prod.xtc in $D"

log "PBC cleanup (whole -> nojump -> centred compact)"
printf 'SOLU\n'        | gmx trjconv -s prod.tpr -f prod.xtc -n index.ndx -pbc whole  -o .w.xtc
printf 'SOLU\n'        | gmx trjconv -s prod.tpr -f .w.xtc   -n index.ndx -pbc nojump -o .nj.xtc
printf 'SOLU\nSOLU\n'  | gmx trjconv -s prod.tpr -f .nj.xtc  -n index.ndx -center \
                                     -pbc mol -ur compact -o "$A/clean.xtc"
printf 'SOLU\n'        | gmx trjconv -s prod.tpr -f "$A/clean.xtc" -n index.ndx \
                                     -dump 0 -o "$A/clean.gro"
rm -f .w.xtc .nj.xtc

cd "$A"
# Protein systems use Backbone / C-alpha; the free hairpin (no protein) uses
# the RNA sugar-phosphate backbone and its phosphorus atoms instead.
if grep -q '\[ Backbone \]' ../index.ndx; then BB=Backbone; CA=C-alpha; else BB=RNA_BB; CA=RNA_P; fi
log "RMSD / Rg / RMSF / SASA / DSSP  (groups: $BB, $CA)"
printf '%s\n%s\n' "$BB" "$BB" | gmx rms   -s ../prod.tpr -f clean.xtc -n ../index.ndx -o rmsd_backbone.xvg -tu ns
printf '%s\n' "$BB"           | gmx gyrate -s ../prod.tpr -f clean.xtc -n ../index.ndx -o gyrate.xvg
# discard the first skip_ns (registry) before RMSF -- fluctuations during
# relaxation are not the fluctuations you want to report
printf '%s\n' "$CA"            | gmx rmsf  -s ../prod.tpr -f clean.xtc -n ../index.ndx -o rmsf_ca.xvg -res -b "$SKIP_PS"
printf 'SOLU\n'               | gmx sasa  -s ../prod.tpr -f clean.xtc -n ../index.ndx -o sasa.xvg
gmx dssp -s ../prod.tpr -f clean.xtc -n ../index.ndx -o dssp.dat -num dssp_num.xvg 2>/dev/null \
  || log "gmx dssp unavailable (needs GROMACS >= 2023); fall back to do_dssp"

log "PCA (essential dynamics), after the first $(( SKIP_PS / 1000 )) ns"
printf '%s\n%s\n' "$CA" "$CA" | gmx covar  -s ../prod.tpr -f clean.xtc -n ../index.ndx \
                                         -o eigenval.xvg -v eigenvec.trr -b "$SKIP_PS" >/dev/null
printf '%s\n%s\n' "$CA" "$CA" | gmx anaeig -s ../prod.tpr -f clean.xtc -n ../index.ndx \
                                         -v eigenvec.trr -first 1 -last 1 -proj pc1.xvg -b "$SKIP_PS" >/dev/null

if [ "$HAS_RNA" = "yes" ]; then
  log "RNA-specific analysis"
  printf 'RNA\nRNA\n' | gmx rms -s ../prod.tpr -f clean.xtc -n ../index.ndx -o rmsd_rna.xvg -tu ns \
    || log "no RNA group in index.ndx -- add one with make_ndx"
  # Hydrogen bonds: protein-RNA when there is a protein, within the RNA when not.
  # GROMACS >= 2024 rewrote `gmx hbond` (selections via -r/-t) and kept the old
  # tool as `hbond-legacy`. Try the new interface, then legacy, then the
  # pre-2024 `hbond`. None of them is fatal: the count is a cross-check only.
  if grep -q '\[ Protein \]' ../index.ndx; then
    HB_SEL=(-r 'group "Protein"' -t 'group "RNA"'); HB_GRP='Protein\nRNA\n'
  else
    HB_SEL=(-r 'group "RNA"');                      HB_GRP='RNA\nRNA\n'
  fi
  HB_IN=(-s ../prod.tpr -f clean.xtc -n ../index.ndx -num hbond_num.xvg)
  if   gmx hbond "${HB_IN[@]}" "${HB_SEL[@]}" < /dev/null > hbond.log 2>&1; then :
  elif printf "$HB_GRP" | gmx hbond-legacy "${HB_IN[@]}" >> hbond.log 2>&1; then :
  elif printf "$HB_GRP" | gmx hbond "${HB_IN[@]}" >> hbond.log 2>&1; then :
  else log "hydrogen-bond count skipped (see $A/hbond.log); everything else is unaffected"
  fi
  # ---- Mg2+ sanity check ------------------------------------------------
  # Mg2+ does not equilibrate on this timescale. If an ion has parked itself
  # in the protein-RNA interface, every energy downstream is contaminated.
  log "Mg2+ interface check -- inspect mindist_mg.xvg by eye"
  printf 'MG\nSOLU\n' | gmx mindist -s ../prod.tpr -f clean.xtc -n ../index.ndx \
                                    -od mindist_mg.xvg > mindist.log 2>&1 \
    || log "Mg2+ distance check skipped (no MG group, or see $A/mindist.log)"
fi

log "convergence"
pyrun "$REPO/analysis/convergence.py" block  rmsd_backbone.xvg | tail -4
pyrun "$REPO/analysis/convergence.py" cosine pc1.xvg
log "done: $A"
