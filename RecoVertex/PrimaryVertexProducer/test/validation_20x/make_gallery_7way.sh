#!/bin/bash
# Seven-way gallery from a CPU job and a GPU job (each from ./run_validation.sh). The GPU job also contains the CPU-only
# collections; only its Alpaka and GNN folders (computed on the GPU) are used here.
#   usage: ./make_gallery_7way.sh <harvested_cpu.root> <harvested_gpu.root> <outdir>
C=$1; G=$2; OUT=${3:-gallery_7way}
[ -f "$C" ] && [ -f "$G" ] || { echo "usage: $0 <harvested_cpu.root> <harvested_gpu.root> <outdir>"; exit 1; }
mkdir -p "$OUT"
python3 "$(dirname "$0")/render_dqm_gallery.py" \
  --old "$C" --old-collection 3D --inblocks "$C" --inblocks-collection InBlocks \
  --alpaka "$C" --alpaka-collection Alpaka --alpaka-cpu "$C" --alpaka-cpu-collection 4D \
  --extra "$C" --extra-collection GNN --extra2 "$G" --extra2-collection Alpaka --extra3 "$G" --extra3-collection GNN \
  --base-dir "DQMData/Run 1/MTD/Run summary/Vertices/{coll}" \
  --display "Old=Legacy D.A. (3D)" --display "InBlocks=D.A. In Blocks (C++)" --display "Alpaka=D.A. In Blocks (Alpaka CPU)" \
  --display "Alpaka-CPU=Legacy D.A. (4D)" --display "Extra=GNN (CPU)" --display "Extra2=D.A. In Blocks (Alpaka GPU)" --display "Extra3=GNN (GPU)" \
  --out "$OUT/pv_dqm_gallery.html" --workdir "$OUT/work" --title "MTD vertex validation: DA in blocks and GNN, CPU and GPU"
