#!/bin/bash
# 4 beam-spot-constrained runs (legacy, C++ in blocks, Alpaka CPU, Alpaka GPU), 91 events ZMM, harvested.
OUT=/home/claude/scratch/gallery_16_bs
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc13 SITECONFIG_PATH=/home/claude/siteconf
cd /home/claude/CMSSW_16_1_0_pre4/src && eval `scram runtime -sh`
cd RecoVertex/PrimaryVertexProducer/test
run() { # tag mode backend
  tag=$1; mode=$2; be=$3
  args="-m $mode --bs -d -i Zmumu_16X_local.txt -o $OUT/$tag.root"
  [ -n "$be" ] && args="$args -b $be"
  cmsRun testVertexing.py $args >& $OUT/$tag.log
  echo "$tag: events=$(grep -c 'Begin processing' $OUT/$tag.log) crash=$(grep -ciE 'out of range|cudaError|Assertion|Fatal Exception|segmentation' $OUT/$tag.log)"
  echo "file:$OUT/${tag}_DQM.root" > $OUT/${tag}_list.txt
  rm -f DQM_V0001_R000000001__Global__CMSSW_*__RECO.root
  cmsRun harvester.py inputFileList=$OUT/${tag}_list.txt outputFile=unused.root >& $OUT/harvest_$tag.log
  mv DQM_V0001_R000000001__Global__CMSSW_*__RECO.root $OUT/${tag}_harv.root
}
run old Old ""
run inblocks InBlocks "" # (rerun with --minndof 2, see the gallery subtitle)
run alpaka_cpu Alpaka cpu
run alpaka_gpu Alpaka gpu-nvidia
echo ALLDONE
