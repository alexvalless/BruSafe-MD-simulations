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
log "RMSD / Rg / RMSF / SASA / DSSP"
printf 'Backbone\nBackbone\n' | gmx rms   -s ../prod.tpr -f clean.xtc -n ../index.ndx -o rmsd_backbone.xvg -tu ns
printf 'Backbone\n'           | gmx gyrate -s ../prod.tpr -f clean.xtc -n ../index.ndx -o gyrate.xvg
# discard the first 50 ns before RMSF -- fluctuations during relaxation are not
# the fluctuations you want to report
printf 'C-alpha\n'            | gmx rmsf  -s ../prod.tpr -f clean.xtc -n ../index.ndx -o rmsf_ca.xvg -res -b 50000
printf 'SOLU\n'               | gmx sasa  -s ../prod.tpr -f clean.xtc -n ../index.ndx -o sasa.xvg
gmx dssp -s ../prod.tpr -f clean.xtc -n ../index.ndx -o dssp.dat -num dssp_num.xvg 2>/dev/null \
  || log "gmx dssp unavailable (needs GROMACS >= 2023); fall back to do_dssp"

log "PCA (essential dynamics), last 200 ns"
printf 'C-alpha\nC-alpha\n' | gmx covar  -s ../prod.tpr -f clean.xtc -n ../index.ndx \
                                         -o eigenval.xvg -v eigenvec.trr -b 50000 >/dev/null
printf 'C-alpha\nC-alpha\n' | gmx anaeig -s ../prod.tpr -f clean.xtc -n ../index.ndx \
                                         -v eigenvec.trr -first 1 -last 1 -proj pc1.xvg -b 50000 >/dev/null

if [ "$HAS_RNA" = "yes" ]; then
  log "RNA-specific analysis"
  printf 'RNA\nRNA\n' | gmx rms -s ../prod.tpr -f clean.xtc -n ../index.ndx -o rmsd_rna.xvg -tu ns \
    || log "no RNA group in index.ndx -- add one with make_ndx"
  printf 'SOLU\n'     | gmx hbond -s ../prod.tpr -f clean.xtc -n ../index.ndx -num hbond_num.xvg \
    || true
  # ---- Mg2+ sanity check ------------------------------------------------
  # Mg2+ does not equilibrate on this timescale. If an ion has parked itself
  # in the protein-RNA interface, every energy downstream is contaminated.
  log "Mg2+ interface check -- inspect mindist_mg.xvg by eye"
  printf 'MG\nSOLU\n' | gmx mindist -s ../prod.tpr -f clean.xtc -n ../index.ndx \
                                    -od mindist_mg.xvg 2>/dev/null \
    || log "no MG group -- skipping (fine for apo systems)"
fi

log "convergence"
python3 "$REPO/analysis/convergence.py" block  rmsd_backbone.xvg | tail -4
python3 "$REPO/analysis/convergence.py" cosine pc1.xvg
log "done: $A"
