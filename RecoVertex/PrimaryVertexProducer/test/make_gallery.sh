#!/bin/bash
# Run Old / InBlocks / Alpaka over the same input sample, harvest DQM, and
# build a single-file HTML comparison gallery.
#
# Usage (run from RecoVertex/PrimaryVertexProducer/test, with cmsenv already
# sourced and the package already built):
#   ./make_gallery.sh [input_file_list] [alpaka_backend] [outdir]
#
#   input_file_list : passed to testVertexing.py -i (default: Zmumu_16X.txt)
#   alpaka_backend   : passed to testVertexing.py -b for the Alpaka run
#                      (default: auto; use gpu-nvidia on a GPU node)
#   outdir           : where to put intermediate ROOT files + the gallery
#                      (default: ./gallery_run)
#
# Any of Old/InBlocks/Alpaka can fail or be skipped (e.g. no GPU available)
# -- the gallery is built from whichever ones succeeded, as long as at
# least two did.

set -uo pipefail

INPUT=${1:-Zmumu_16X.txt}
BACKEND=${2:-auto}
OUTDIR=${3:-gallery_run}

mkdir -p "$OUTDIR"
cd "$(dirname "$0")"

check_proxy() {
  echo "=== Checking grid proxy ==="
  if voms-proxy-info --exists --valid 1:0 >/dev/null 2>&1; then
    echo "   OK ($(voms-proxy-info --timeleft 2>/dev/null)s remaining)"
    return 0
  fi
  echo "   No valid grid proxy (or expiring within an hour) -- running voms-proxy-init."
  echo "   This needs your grid certificate passphrase, entered interactively below."
  voms-proxy-init --rfc --voms cms -valid 192:00
  if [ $? -ne 0 ]; then
    echo "!! voms-proxy-init failed -- can't read input files over xrootd without it." >&2
    exit 1
  fi
}

# Run a command in the background, streaming its log live via `tail -f`
# (which exits on its own once the job's PID is gone -- GNU coreutils
# --pid), then recover the job's real exit code via `wait`.
run_and_tail() {
  local log=$1
  shift
  "$@" >& "$log" &
  local pid=$!
  tail -f --pid="$pid" "$log" 2>/dev/null
  wait "$pid"
  return $?
}

run_mode() {
  local mode=$1
  local extra_args=$2
  local tag=$3
  echo "=== Running $mode ($tag) ==="
  run_and_tail "$OUTDIR/${tag}.log" \
    cmsRun testVertexing.py -m "$mode" -d -i "$INPUT" $extra_args -o "$OUTDIR/${tag}.root"
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "!! $mode run failed, see $OUTDIR/${tag}.log" >&2
    return 1
  fi
  if grep -q "segmentation violation" "$OUTDIR/${tag}.log"; then
    echo "!! $mode run crashed (segfault), see $OUTDIR/${tag}.log" >&2
    return 1
  fi
  echo "   -> $OUTDIR/${tag}_DQM.root"
  return 0
}

harvest() {
  local tag=$1
  echo "=== Harvesting $tag ==="
  echo "file:$OUTDIR/${tag}_DQM.root" > "$OUTDIR/${tag}_DQM_list.txt"
  rm -f DQM_V0001_R000000001__Global__CMSSW_*__RECO.root
  run_and_tail "$OUTDIR/harvest_${tag}.log" \
    cmsRun harvester.py inputFileList="$OUTDIR/${tag}_DQM_list.txt" outputFile=unused.root
  local rc=$?
  if [ $rc -ne 0 ] || ! ls DQM_V0001_R000000001__Global__CMSSW_*__RECO.root >/dev/null 2>&1; then
    echo "!! Harvest of $tag failed, see $OUTDIR/harvest_${tag}.log" >&2
    return 1
  fi
  mv DQM_V0001_R000000001__Global__CMSSW_*__RECO.root "$OUTDIR/${tag}_harvested.root"
  echo "   -> $OUTDIR/${tag}_harvested.root"
  return 0
}

check_proxy

HAVE_OLD=0
HAVE_INBLOCKS=0
HAVE_ALPAKA=0

if run_mode Old "" old; then
  harvest old && HAVE_OLD=1
fi

if run_mode InBlocks "" inblocks; then
  harvest inblocks && HAVE_INBLOCKS=1
fi

if run_mode Alpaka "-b $BACKEND" alpaka; then
  harvest alpaka && HAVE_ALPAKA=1
fi

ARGS=()
[ "$HAVE_OLD" = "1" ] && ARGS+=(--old "$OUTDIR/old_harvested.root")
[ "$HAVE_INBLOCKS" = "1" ] && ARGS+=(--inblocks "$OUTDIR/inblocks_harvested.root")
[ "$HAVE_ALPAKA" = "1" ] && ARGS+=(--alpaka "$OUTDIR/alpaka_harvested.root")

NSUCCESS=$((HAVE_OLD + HAVE_INBLOCKS + HAVE_ALPAKA))
if [ "$NSUCCESS" -lt 2 ]; then
  echo "!! Fewer than two modes succeeded ($NSUCCESS) -- can't build a comparison gallery." >&2
  exit 1
fi

echo "=== Building gallery from $NSUCCESS mode(s) ==="
python3 render_dqm_gallery.py "${ARGS[@]}" \
  --out "$OUTDIR/pv_dqm_gallery.html" \
  --workdir "$OUTDIR/gallery_pngs" \
  --subtitle "input: $INPUT, Alpaka backend: $BACKEND"

echo "=== Done: $OUTDIR/pv_dqm_gallery.html ==="
