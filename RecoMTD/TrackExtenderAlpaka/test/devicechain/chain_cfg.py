# MTD timing + Alpaka 3D vertices (+ TOF PID) on the skim of skim_cfg.py; writes the products for exact comparisons.
# usage: cmsRun chain_cfg.py inputFiles=file:skim.root mode=new|tracks|dev2 output=file.root [maxEvents=N] [threads=N]
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Phase2C26I13M9_cff import Phase2C26I13M9

options = VarParsing('analysis')
options.register('mode', 'new', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'new, tracks or dev2')
options.register('output', 'chain.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'output')
options.register('threads', 4, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'threads')
options.register('backend', 'cpu', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'cpu or gpu-nvidia')
options.parseArguments()

process = cms.Process('CHAIN', Phase2C26I13M9)
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D128Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Accelerators_cff')
process.load('HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

process.source = cms.Source('PoolSource', fileNames=cms.untracked.vstring(options.inputFiles))
process.maxEvents.input = options.maxEvents
process.options.accelerators = [options.backend]
process.options.numberOfThreads = options.threads
process.options.wantSummary = True
process.MessageLogger.cerr.FwkReport.reportEvery = 10

# MTD rechits rebuilt from the stored clusters (same CPE as the reconstruction)
process.load('RecoLocalFastTime.FTLClusterizer.MTDCPEESProducers_cff')
process.load('RecoLocalFastTime.FTLRecProducers.MTDTimeCalibESProducers_cff')
from RecoLocalFastTime.FTLRecProducers.mtdTrackingRecHits_cfi import mtdTrackingRecHits as _mtdTRH
process.mtdTrackingRecHits = _mtdTRH.clone()
hits = cms.InputTag('mtdTrackingRecHits', '', 'CHAIN')

process.trackExtenderWithMTD = cms.EDProducer('MtdTimingToValueMaps',
    timingSrc=cms.InputTag('mtdTimingProducer'), trackClones=cms.bool(True), hitsSrc=hits)
from RecoMTD.TimingIDTools.mtdTrackQualityMVA_cfi import mtdTrackQualityMVA as _mva
process.mtdTrackQualityMVA = _mva.clone()

import RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesAlpaka_cfi as pv
for name in ['pvVertexSoA', 'pvVertexSoAWithBS', 'unsortedOfflinePrimaryVerticesFromSoA',
             'unsortedOfflinePrimaryVerticesFromSoAWithBS']:
    setattr(process, name, getattr(pv, name).clone())

if options.mode in ('new', 'tracks'):
    # one packed input for the whole device chain, BTL geometry from the event setup
    process.mtdBtlGeometryESProducer = cms.ESProducer('MtdBtlGeometryESProducer@alpaka')
    process.mtdEventInput = cms.EDProducer('MtdEventInputProducer@alpaka', hitsSrc=hits,
        TkFilterParameters=pv.pvTracksSoA.TkFilterParameters.clone())
    process.mtdTimingProducer = cms.EDProducer('MtdTimingProducer@alpaka', inputSrc=cms.InputTag('mtdEventInput'))
    tasks = [process.mtdEventInput, process.mtdTimingProducer]
    if options.mode == 'new':
        for name in ['pvVertexSoA', 'pvVertexSoAWithBS']:
            getattr(process, name).TrackLabel = 'mtdEventInput'
            getattr(process, name).inputType = 'mtdEventInput'
    else:
        # vertices from the separate track SoA of PortableTrackSoAProducer (regression of that path)
        process.pvTracksSoA = pv.pvTracksSoA.clone()
        tasks.append(process.pvTracksSoA)
elif options.mode == 'dev2':
    # step 2: one clustering with both fits, vertex time and TOF PID on the device; the host versions run
    # next to them for the comparison (host BDT from the maps, host TOF PID on the converted vertices)
    process.mtdBtlGeometryESProducer = cms.ESProducer('MtdBtlGeometryESProducer@alpaka')
    process.mtdEventInput = cms.EDProducer('MtdEventInputProducer@alpaka', hitsSrc=hits,
        TkFilterParameters=pv.pvTracksSoA.TkFilterParameters.clone())
    process.mtdTimingProducer = cms.EDProducer('MtdTimingProducer@alpaka', inputSrc=cms.InputTag('mtdEventInput'))
    process.pvVertexSoA.TrackLabel = 'mtdEventInput'
    process.pvVertexSoA.inputType = 'mtdEventInput'
    process.pvVertexSoA.alsoWithBS = cms.bool(True)
    process.pvVertexSoA.timingSrc = cms.InputTag('mtdTimingProducer')
    del process.pvVertexSoAWithBS
    process.unsortedOfflinePrimaryVerticesFromSoA.vertexTimeParameters.algorithm = 'fromSoA'
    process.unsortedOfflinePrimaryVerticesFromSoAWithBS.soaVertex = cms.InputTag('pvVertexSoA', 'WithBS')
    process.unsortedOfflinePrimaryVerticesFromSoAWithBS.vertexTimeParameters.algorithm = 'fromSoA'
    process.tofPidProducer = cms.EDProducer('TofPidProducer@alpaka', inputSrc=cms.InputTag('mtdEventInput'),
        timingSrc=cms.InputTag('mtdTimingProducer'), vtxsSrc=cms.InputTag('pvVertexSoA'))
    process.tofPID4DnoPIDDevice = cms.EDProducer('TofPidToValueMaps', tofPidSrc=cms.InputTag('tofPidProducer'))
    from RecoMTD.TimingIDTools.tofPIDProducer_cfi import tofPIDProducer as _tofPID
    process.tofPID4DnoPID = _tofPID.clone(vtxsSrc='unsortedOfflinePrimaryVerticesFromSoA', MVASel=False,
        probPion=1., probKaon=1., probProton=1., vertexReassignment=True, vtxMaxSigmaT=0.025, maxDz=0.1,
        maxDtSignificance=5., minProbHeavy=0.75, fixedT0Error=0.)
    tasks = [process.mtdEventInput, process.mtdTimingProducer, process.tofPidProducer,
             process.tofPID4DnoPIDDevice, process.tofPID4DnoPID]
else:
    raise ValueError('mode must be new, tracks or dev2')

process.chainTask = cms.Task(process.mtdTrackingRecHits, process.trackExtenderWithMTD, process.mtdTrackQualityMVA,
                             process.pvVertexSoA, process.unsortedOfflinePrimaryVerticesFromSoA,
                             process.unsortedOfflinePrimaryVerticesFromSoAWithBS, *tasks)
if hasattr(process, 'pvVertexSoAWithBS'):
    process.chainTask.add(process.pvVertexSoAWithBS)
process.out = cms.OutputModule('PoolOutputModule', fileName=cms.untracked.string(options.output),
    outputCommands=cms.untracked.vstring('drop *',
        'keep *_trackExtenderWithMTD_*_CHAIN', 'keep *_mtdTrackQualityMVA_*_CHAIN',
        'keep recoVertexs_unsortedOfflinePrimaryVerticesFromSoA*_*_CHAIN', 'keep *_tofPID4DnoPID*_*_CHAIN'))
process.e = cms.EndPath(process.out, process.chainTask)
