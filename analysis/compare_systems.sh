#!/usr/bin/env bash
# Cross-system comparisons.
#
#   ./compare_systems.sh <system_A> <system_B>
#
# The headline use: compare_systems.sh S1_wt_cc_apo S2_sccp_apo
# Subspace overlap between the two essential subspaces is the rigorous way to
# claim "the linker does not alter the dimer's dynamics". 
set -euo pipefail
source "$(dirname "$0")/../scripts/lib.sh"
need gmx

A="${1:?usage: compare_systems.sh <system_A> <system_B>}"
B="${2:?usage: compare_systems.sh <system_A> <system_B>}"
OUT="$REPO/analysis_out/${A}__vs__${B}"
mkdir -p "$OUT"; cd "$OUT"

log "essential subspace overlap (rep1 of each)"
gmx anaeig -v "$RUNS/$A/rep1/analysis/eigenvec.trr" \
           -v2 "$RUNS/$B/rep1/analysis/eigenvec.trr" \
           -s  "$RUNS/$A/rep1/prod.tpr" \
           -first 1 -last 10 -over overlap.xvg
log "overlap.xvg written -- >0.7 over the first 10 eigenvectors supports"
log "'same essential dynamics'; <0.5 means the linker changed something real"

log "replica agreement within each system"
python3 "$REPO/analysis/convergence.py" replicas "$RUNS/$A"/rep*/analysis/rmsf_ca.xvg \
  > rmsf_agreement_"$A".txt
python3 "$REPO/analysis/convergence.py" replicas "$RUNS/$B"/rep*/analysis/rmsf_ca.xvg \
  > rmsf_agreement_"$B".txt
cat rmsf_agreement_"$A".txt rmsf_agreement_"$B".txt

log "IMPORTANT: the difference between A and B is only meaningful if it"
log "exceeds the inter-replica spread within each. Check that before claiming"
log "the linker does anything."
