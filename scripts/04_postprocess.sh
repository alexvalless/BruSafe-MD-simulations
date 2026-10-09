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
  if   gmx_timeout 300 hbond "${HB_IN[@]}" "${HB_SEL[@]}" < /dev/null > hbond.log 2>&1; then :
  elif printf "$HB_GRP" | gmx_timeout 300 hbond-legacy "${HB_IN[@]}" >> hbond.log 2>&1; then :
  elif printf "$HB_GRP" | gmx_timeout 300 hbond "${HB_IN[@]}" >> hbond.log 2>&1; then :
  else log "hydrogen-bond count skipped (see $A/hbond.log); everything else is unaffected"
  fi
  # ---- Mg2+ sanity check ------------------------------------------------
  # Mg2+ does not equilibrate on this timescale. If an ion has parked itself
  # in the protein-RNA interface, every energy downstream is contaminated.
  # Uses the RAW prod.xtc: clean.xtc holds only the SOLU group, and the ions
  # (MG) are not in it, so any MG group is beyond the end of that trajectory
  # (the old `gmx mindist -f clean.xtc` read out of range and segfaulted).
  # `gmx pairdist` takes the minimum image for every pair, so unwrapped
  # molecules are fine. It gives the minimum distance of EACH ion to the RNA
  # and to the protein.
  log "Mg2+ interface check -- per-ion distance to RNA and to protein"
  if grep -q '\[ MG \]' ../index.ndx; then
    : > mindist.log
    gmx_timeout 900 pairdist -s ../prod.tpr -f ../prod.xtc -n ../index.ndx \
        -ref 'group "MG"' -sel 'group "RNA"' -refgrouping none -type min \
        -o mindist_mg_rna.xvg >> mindist.log 2>&1 \
      || log "Mg2+ to RNA distances skipped (see $A/mindist.log)"
    if grep -q '\[ Protein \]' ../index.ndx; then
      gmx_timeout 900 pairdist -s ../prod.tpr -f ../prod.xtc -n ../index.ndx \
          -ref 'group "MG"' -sel 'group "Protein"' -refgrouping none -type min \
          -o mindist_mg_protein.xvg >> mindist.log 2>&1 \
        || log "Mg2+ to protein distances skipped (see $A/mindist.log)"
      if [ -s mindist_mg_rna.xvg ] && [ -s mindist_mg_protein.xvg ]; then
        # columns: time, then one distance (nm) per ion; both files have the same frames
        paste <(grep -v '^[#@]' mindist_mg_rna.xvg | tr -d '\r') \
              <(grep -v '^[#@]' mindist_mg_protein.xvg | tr -d '\r') \
          | awk -v cut=0.5 '
              { n = NF / 2 - 1; N = n; F = NR
                for (i = 1; i <= n; i++) {
                  r = $(1 + i); p = $(n + 2 + i)
                  if (r < cut) nr[i]++; if (p < cut) np[i]++; if (r < cut && p < cut) both[i]++ } }
              END { printf "# fraction of %d frames with the ion within %.1f nm of the RNA / the protein / BOTH\n", F, cut
                    print "# ion    RNA  protein   BOTH   (BOTH > 0: the ion sits at the protein-RNA interface)"
                    for (i = 1; i <= N; i++) printf "%5d  %5.2f  %7.2f  %5.2f\n", i, nr[i]/F, np[i]/F, both[i]/F }' \
          > mg_interface.txt
        log "Mg2+ interface table: $A/mg_interface.txt"
      fi
    fi
  else
    log "no MG group -- skipping (fine for apo systems)"
  fi
fi

log "convergence"
pyrun "$REPO/analysis/convergence.py" block  rmsd_backbone.xvg | tail -4
pyrun "$REPO/analysis/convergence.py" cosine pc1.xvg
log "done: $A"
