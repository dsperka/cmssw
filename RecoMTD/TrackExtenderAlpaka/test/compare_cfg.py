# Run the new MTD timing (Alpaka) on the skim and dump legacy vs new vs truth per track.
# usage: cmsRun compare_cfg.py [inputFile] [outputCsv] [backend: cpu|gpu-nvidia]
import sys
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C26I13M9_cff import Phase2C26I13M9

args = [a for a in sys.argv[1:] if not a.endswith('.py')]
inFile = args[0] if len(args) > 0 else 'file:/home/claude/scratch/mtd_alpaka_dev/skim/mtd_skim.root'
outCsv = args[1] if len(args) > 1 else 'compare.csv'
backend = args[2] if len(args) > 2 else 'cpu'

process = cms.Process('MTDNEW', Phase2C26I13M9)
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D128Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Accelerators_cff')
process.load('HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

process.source = cms.Source('PoolSource', fileNames=cms.untracked.vstring(inFile))
process.maxEvents.input = int(__import__("os").environ.get("MTDT_NEV","-1"))
process.options.accelerators = [backend]
process.options.wantSummary = True
process.MessageLogger.cerr.FwkReport.reportEvery = 10

# MTD tracking rechits do not keep their local position/error in the file: rebuild them from the
# stored clusters, exactly as the reconstruction does (same CPE).
process.load('RecoLocalFastTime.FTLClusterizer.MTDCPEESProducers_cff')
process.load('RecoLocalFastTime.FTLRecProducers.MTDTimeCalibESProducers_cff')
from RecoLocalFastTime.FTLRecProducers.mtdTrackingRecHits_cfi import mtdTrackingRecHits as _mtdTRH
process.mtdTrackingRecHits = _mtdTRH.clone()
process.mtdTimingProducer = cms.EDProducer('MtdTimingProducer@alpaka',
    hitsSrc=cms.InputTag('mtdTrackingRecHits', '', 'MTDNEW'))
process.trackExtenderWithMTDAlpaka = cms.EDProducer('MtdTimingToValueMaps',
    timingSrc=cms.InputTag('mtdTimingProducer'), trackClones=cms.bool(True),
    hitsSrc=cms.InputTag('mtdTrackingRecHits', '', 'MTDNEW'))
# truth: joined offline from the earlier hit-based truth dump (see analysis/step5_check.py)
process.compare = cms.EDAnalyzer('MtdTimingCompare',
    tracks=cms.InputTag('generalTracks'),
    recoToSim=cms.InputTag('none'),
    legacy=cms.string('trackExtenderWithMTD'),
    alpaka=cms.string('trackExtenderWithMTDAlpaka'),
    fileName=cms.string(outCsv))
process.p = cms.Path(process.mtdTrackingRecHits + process.mtdTimingProducer + process.trackExtenderWithMTDAlpaka + process.compare)
