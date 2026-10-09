#!/usr/bin/env bash
# Pack the FINISHED replicas of THIS machine into archives that are easy to
# upload (OneDrive, Zenodo) and to verify. Safe to rerun: replicas that are
# already packed are skipped.
#
#   ./scripts/pack_results.sh            # pack, and copy to OneDrive if found
#   ./scripts/pack_results.sh list       # show what would be packed, do nothing
#
# Per finished replica (prod.gro exists) it writes, into pack/<hostname>/:
#   <system>_rep<N>_light.tar.gz   everything needed to reproduce and analyse:
#                                  topology, index, mdp, seed, PROVENANCE,
#                                  prod.tpr/.gro/.edr/.log and the whole
#                                  analysis/ folder (incl. the dry clean.xtc).
#                                  Small (tens of MB): this is the one for Zenodo.
#   <system>_rep<N>_traj.tar       prod.xtc (full system with water) + prod.cpt:
#                                  big (about 1 GB), for backup only.
#   SHA256SUMS.txt                 checksums of every archive.
#
# If a OneDrive folder is found (~/OneDrive*), the archives are copied to
# <OneDrive>/BruSafe-MD/<hostname>/ ; otherwise upload the pack/<hostname>
# folder from the browser (onedrive.live.com / the Tec Microsoft 365 page).
# Override the destination with BRUSAFE_PACK_DEST=/path.
set -euo pipefail
source "$(dirname "$0")/lib.sh"

MODE="${1:-pack}"
HOST="$(hostname)"
OUTDIR="$REPO/pack/$HOST"
mkdir -p "$OUTDIR"

DEST="${BRUSAFE_PACK_DEST:-}"
if [ -z "$DEST" ]; then
  for od in "$HOME"/OneDrive*; do
    [ -d "$od" ] && { DEST="$od/BruSafe-MD/$HOST"; break; }
  done
fi

n_done=0; n_new=0
for d in "$RUNS"/*/rep*; do
  [ -d "$d" ] && [ -f "$d/prod.gro" ] || continue
  sys="$(basename "$(dirname "$d")")"; rep="${d##*rep}"
  light="$OUTDIR/${sys}_rep${rep}_light.tar.gz"
  traj="$OUTDIR/${sys}_rep${rep}_traj.tar"
  n_done=$((n_done + 1))
  if [ -f "$light" ] && [ -f "$traj" ] && [ "$light" -nt "$d/prod.gro" ]; then
    echo "  already packed: $sys rep$rep"; continue
  fi
  if [ "$MODE" = list ]; then echo "  would pack:     $sys rep$rep"; continue; fi

  log "packing $sys rep$rep"
  # light: everything except the heavy trajectory/checkpoints and GROMACS backups
  tar -czf "$light" -C "$(dirname "$d")" \
      --exclude='#*#' --exclude='prod.xtc' --exclude='prod.cpt' --exclude='prod_prev.cpt' \
      --exclude='*.trr' "rep$rep"
  # traj: the big files
  files=()
  for f in prod.xtc prod.cpt; do [ -f "$d/$f" ] && files+=("rep$rep/$f"); done
  tar -cf "$traj" -C "$(dirname "$d")" "${files[@]}"
  n_new=$((n_new + 1))
done

if [ "$MODE" = list ]; then echo "$n_done finished replica(s) on this machine"; exit 0; fi

if [ "$n_done" -eq 0 ]; then log "no finished replica on this machine yet"; exit 0; fi

( cd "$OUTDIR" && sha256sum ./*.tar* | sed 's| \./| |' | tr -d '\r' > SHA256SUMS.txt )
log "$n_done finished replica(s), $n_new newly packed -> $OUTDIR"
du -sh "$OUTDIR"/* | sed 's/^/  /'

if [ -n "$DEST" ]; then
  mkdir -p "$DEST"
  # small files first: if the connection drops, the important ones are in
  cp -u "$OUTDIR"/*_light.tar.gz "$OUTDIR"/SHA256SUMS.txt "$DEST"/
  cp -u "$OUTDIR"/*_traj.tar "$DEST"/
  log "copied to $DEST (OneDrive will sync it; check the icon says 'up to date')"
else
  log "no OneDrive folder found: upload $OUTDIR from the browser"
fi
