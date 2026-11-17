# Skim for chain_cfg.py: the full reconstruction of ../validation/mtdValidation_cfg.py with the legacy MTD
# timing, keeping only what the device chain and its comparisons read (tracks, MTD clusters and rechits,
# beam spot, the legacy MTD products and the 4D vertices).
#   cmsRun skim_cfg.py inputFiles=file:/path/to/GEN-SIM-RECO.root maxEvents=100 output=skim.root [threads=8]
import os
import sys

import FWCore.ParameterSet.Config as cms

sys.argv.append('timing=legacy')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'validation'))
from mtdValidation_cfg import process, options  # noqa: E402

for name in ['validation_step', 'DQMoutput_step']:
    process.schedule.remove(getattr(process, name))
    delattr(process, name)
del process.DQMoutput

process.skimOutput = cms.OutputModule('PoolOutputModule',
    fileName=cms.untracked.string('file:' + options.output),
    outputCommands=cms.untracked.vstring(
        'drop *',
        'keep *_generalTracks_*_*',
        'keep *_mtdTrackingRecHits_*_*',
        'keep *_mtdClusters_*_*',
        'keep *_offlineBeamSpot_*_*',
        'keep *_trackExtenderWithMTD_*_*',
        'keep *_mix_MergedTrackTruth_*',
        'keep *_trackingParticleRecoTrackAsssociation_*_*',
        'keep *_offlinePrimaryVertices4D_*_*'))
process.skim_step = cms.EndPath(process.skimOutput, process.mtdValidTruthTask)
process.schedule.append(process.skim_step)
