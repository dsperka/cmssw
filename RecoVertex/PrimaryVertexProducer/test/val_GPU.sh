#!/bin/bash

eval `scramv1 runtime -sh`

cd ../../../
clear
cmsenv
#scram b clean
#scram b -j 4
cd RecoVertex/PrimaryVertexProducer/test
cmsenv

#cmsRun vertexTest.py n=1000 verbose=False gpu=False > cpu1000.log
#cmsRun vertexTest.py n=1000 verbose=False gpu=True > gpu1000.log

cmsRun vertexTest.py n=1 verbose=True gpu=False > cpu.log
cmsRun vertexTest.py n=1 verbose=True gpu=True > gpu.log

#harvestTrackValidationPlots.py test_dqm_cpu.root -o cpu.root
#harvestTrackValidationPlots.py test_dqm_gpu.root -o gpu.root
#rm -r plots/
#makeTrackValidationPlots.py gpu.root cpu.root --png --extended
