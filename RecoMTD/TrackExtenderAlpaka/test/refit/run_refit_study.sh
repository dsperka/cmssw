#!/bin/bash
# Refit study: legacy (refit with the MTD hit) vs Alpaka (no refit) vs truth, on the 100-event skim of
# ../devicechain. Run inside cmsenv from a work directory:
#
#   run_refit_study.sh skim.root [threads]
#
# 1. the Alpaka MTD timing on the skim (chain_cfg.py mode=new) -> chain.root
# 2. per-track table (make_table.py, FWLite)                    -> table.npz
# 3. tables: refit vs original vs truth momentum, Alpaka vs legacy TOF inputs (analyze_refit.py)
#            t0 vs the true production time (analyze_truth.py)  -> refit.txt, truth.txt
# 4. plots (make_plots.py)                                      -> plots/*.png
set -e
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SKIM=$(readlink -f "$1")
THREADS=${2:-8}
[ -s chain.root ] || cmsRun "$HERE/../devicechain/chain_cfg.py" inputFiles=file:$SKIM mode=new output=chain.root \
  threads=$THREADS > chain.log 2>&1
[ -s table.npz ] || python3 "$HERE/make_table.py" "$SKIM" chain.root table.npz
python3 "$HERE/analyze_refit.py" table.npz | tee refit.txt
python3 "$HERE/analyze_truth.py" table.npz | tee truth.txt
python3 "$HERE/make_plots.py" table.npz plots
