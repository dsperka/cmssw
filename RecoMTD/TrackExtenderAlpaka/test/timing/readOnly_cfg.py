# Read-only job: reads from the input files exactly the branches that the reconstruction job reads
# (the list comes from its framework job report, ReadBranches), and does nothing else.
# Its throughput is the most the disk can give to the reconstruction job.
#   cmsRun readOnly_cfg.py inputFiles=... maxEvents=N branches=branches.txt
import os

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('branches', 'branches.txt', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 'branch names to read, one per line (type_label_instance_process)')
options.parseArguments()

process = cms.Process('READ')
process.source = cms.Source('PoolSource', fileNames=cms.untracked.vstring(options.inputFiles),
                            inputCommands=cms.untracked.vstring('keep *', 'drop *_*_*_RECO'))
process.maxEvents.input = options.maxEvents
process.options.numberOfThreads = 1
process.options.numberOfStreams = 1
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.ThroughputService = cms.untracked.PSet(limit=cms.untracked.int32(10000000))

branches = [b.strip().rstrip('.') for b in open(options.branches) if b.strip()]
process.readAll = cms.EDAnalyzer('GenericConsumer', eventProducts=cms.untracked.vstring(branches))
process.read_step = cms.EndPath(process.readAll)

process.ThroughputService = cms.Service('ThroughputService',
    eventRange=cms.untracked.uint32(100000),
    eventResolution=cms.untracked.uint32(1),
    eventSkip=cms.untracked.uint32(int(os.environ.get('MTD_TIMING_SKIP', '5'))),
    printEventSummary=cms.untracked.bool(True),
    enableDQM=cms.untracked.bool(False))
