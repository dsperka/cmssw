#!/bin/bash
# Reconstruction (all vertex collections) -> MTD validation -> harvesting, on one GEN-SIM-RECO file.
#   usage: ./run_validation.sh <input> [nevents=100] [accelerators=cpu]
#   <input>: a local path, a /store/... LFN, or a root:// URL. accelerators: cpu (default) or gpu-nvidia.
# Run from a scram area where this branch is built (see README.md). Output goes to the current directory.
set -e
IN=$1; N=${2:-100}; ACC=${3:-cpu}
[ -n "$IN" ] || { echo "usage: $0 <input> [nevents] [accelerators]"; exit 1; }
case "$IN" in
  file:*|/store/*|root://*|http*) URL="$IN" ;;
  *) URL="file:$(readlink -f "$IN")" ;;
esac
HERE=$(cd "$(dirname "$0")" && pwd)
sed -e "s#INPUTFILE#${URL}#; s#NEVENTS#${N}#; s#process.options.accelerators = \[\"cpu\"\]#process.options.accelerators = [\"${ACC}\"]#" "$HERE/reco_cfg_template.py" > reco_run.py
cmsRun reco_run.py > reco.log 2>&1
grep "TrigReport Events total" reco.log | head -1
cmsRun "$HERE/mtdValidation_6way_cfg.py" > val.log 2>&1
cmsRun "$HERE/mtdHarvesting_6way_cfg.py" > harv.log 2>&1
ls DQM_V*.root
