#!/bin/bash

harvestTrackValidationPlots.py test_dqm_gpu.root -o gpu.root
rm -r plots/
makeTrackValidationPlots.py gpu.root cpu.root --png --extended
