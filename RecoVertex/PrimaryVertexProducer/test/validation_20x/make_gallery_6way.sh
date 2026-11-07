#!/bin/bash
# Render the 6-way (five collections) gallery from the harvested validation file. Needs ROOT (cmsenv) and render_dqm_gallery.py's dependencies.
#   usage: ./make_gallery_6way.sh <harvested.root> <outdir>
R=${1:-DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root}; G=${2:-gallery_6way}
mkdir -p "$G"
python3 "$(dirname "$0")/render_dqm_gallery.py" \
  --old "$R" --old-collection 3D --inblocks "$R" --inblocks-collection InBlocks \
  --alpaka "$R" --alpaka-collection Alpaka --alpaka-cpu "$R" --alpaka-cpu-collection 4D --extra "$R" --extra-collection GNN \
  --base-dir "DQMData/Run 1/MTD/Run summary/Vertices/{coll}" \
  --display "Old=Legacy D.A. (3D)" --display "InBlocks=D.A. In Blocks (C++)" --display "Alpaka=D.A. In Blocks (Alpaka CPU)" \
  --display "Alpaka-CPU=Legacy D.A. (4D)" --display "Extra=GNN vertex slots" \
  --time-old 775 --time-inblocks 229 --time-alpaka 312 --time-alpaka-cpu 833 --time-extra 1336 \
  --out "$G/pv_dqm_gallery.html" --workdir "$G/work" --title "MTD vertex validation: DA in blocks and GNN on 20_1"
