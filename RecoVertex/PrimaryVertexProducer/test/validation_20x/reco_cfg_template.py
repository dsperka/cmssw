import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C26I13M9_cff import Phase2C26I13M9
from Configuration.ProcessModifiers.vertexSlotGNN_cff import vertexSlotGNN

process = cms.Process('REVTX', Phase2C26I13M9, vertexSlotGNN)

process.load('Configuration.StandardSequences.Services_cff')
process.load('HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi')
process.load('HeterogeneousCore.CUDACore.ProcessAcceleratorCUDA_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D128Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(NEVENTS))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('INPUTFILE'),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop *_SimClusterToCaloParticleAssociation_*_*',
        'drop *_allTrackstersToSimTrackstersAssociationsByLCs_*_*',
        'drop *_allTrackstersToSimTrackstersAssociationsByHits_*_*',
    ),
)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(1),
    wantSummary = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.exe = cms.Path(process.firstStepPrimaryVerticesUnsorted * process.vertexreco)

# ---------------- ported DA-in-blocks developments: InBlocks (C++) and Alpaka ----------------
process.options.accelerators = ["cpu"]

# InBlocks (C++): the standard 3D producer with runInBlocks, fitter as in the 16_X campaign
process.unsortedOfflinePrimaryVerticesInBlocks = process.unsortedOfflinePrimaryVertices.clone(
    TkClusParameters = dict(TkDAClusParameters = dict(runInBlocks = True, block_size = 512, overlap_frac = 0.5, mintrkweight = cms.double(0.25))),
    vertexCollections = cms.VPSet(
        cms.PSet(label = cms.string(""), algorithm = cms.string("WeightedMeanFitter"), chi2cutoff = cms.double(2.5),
                 minNdof = cms.double(0.0), useBeamConstraint = cms.bool(False), maxDistanceToBeam = cms.double(1.0)),
        cms.PSet(label = cms.string("WithBS"), algorithm = cms.string("WeightedMeanFitter"), chi2cutoff = cms.double(2.5),
                 minNdof = cms.double(0.0), useBeamConstraint = cms.bool(True), maxDistanceToBeam = cms.double(1.0))))

# Alpaka: portable tracks -> DA in blocks on the device -> reco::Vertex
process.tracksSoA = cms.EDProducer("PortableTrackSoAProducer@alpaka",
    TrackLabel = cms.InputTag("generalTracks"),
    BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
    TkFilterParameters = cms.PSet(
        maxNormalizedChi2 = cms.double(10.0), minPixelLayersWithHits = cms.int32(2), minSiliconLayersWithHits = cms.int32(5),
        maxD0Significance = cms.double(4.0), maxD0Error = cms.double(1.0), maxDzError = cms.double(1.0), minPt = cms.double(0.0),
        maxEta = cms.double(2.4), trackQuality = cms.string("any"), vertexSize = cms.double(0.006), d0CutOff = cms.double(3.)))
process.beamSpotSoA = cms.EDProducer("BeamSpotDeviceProducer@alpaka", src = cms.InputTag("offlineBeamSpot"))
process.vertexSoA = cms.EDProducer("PrimaryVertexProducerPortable@alpaka",
    TrackLabel = cms.InputTag("tracksSoA"), BeamSpotLabel = cms.InputTag("beamSpotSoA"),
    blockOverlap = cms.double(0.50), blockSize = cms.int32(512),
    TkFitterParameters = cms.PSet(useBeamSpotConstraint = cms.bool(False)),
    TkClusParameters = cms.PSet(coolingFactor = cms.double(0.6), zrange = cms.double(4.), delta_highT = cms.double(1.e-2),
        delta_lowT = cms.double(1.e-3), convergence_mode = cms.int32(0), Tmin = cms.double(2.0), Tpurge = cms.double(2.0),
        Tstop = cms.double(0.5), vertexSize = cms.double(0.006), d0CutOff = cms.double(3.), dzCutOff = cms.double(3.),
        zmerge = cms.double(1e-2), uniquetrkweight = cms.double(0.8), uniquetrkminp = cms.double(0.0), mintrkweight = cms.double(0.25)))
process.unsortedOfflinePrimaryVerticesAlpaka = cms.EDProducer("RecoVertexFromSoAProducer",
    soaVertex = cms.InputTag("vertexSoA"), srcTrack = cms.InputTag("generalTracks"))

def sortedChain(proc, tag):
    # same sorting and timing-PID chain as the standard 3D collection
    setattr(proc, "trackWithVertexRefSelectorBeforeSorting"+tag, proc.trackWithVertexRefSelector.clone(
        vertexTag = "unsortedOfflinePrimaryVertices"+tag, ptMax = 9e99, ptErrorCut = 9e99))
    setattr(proc, "trackRefsForJetsBeforeSorting"+tag, proc.trackRefsForJets.clone(src = "trackWithVertexRefSelectorBeforeSorting"+tag))
    setattr(proc, "offlinePrimaryVertices"+tag, proc.offlinePrimaryVertices.clone(
        vertices = "unsortedOfflinePrimaryVertices"+tag, particles = "trackRefsForJetsBeforeSorting"+tag))
    setattr(proc, "tofPID"+tag, proc.tofPID3D.clone(vtxsSrc = "unsortedOfflinePrimaryVertices"+tag))
    return [getattr(proc, n+tag) for n in ("trackWithVertexRefSelectorBeforeSorting","trackRefsForJetsBeforeSorting","offlinePrimaryVertices","tofPID")]
chainIB = sortedChain(process, "InBlocks")
chainAlpaka = sortedChain(process, "Alpaka")
process.exe2 = cms.Path(process.unsortedOfflinePrimaryVerticesInBlocks + sum(chainIB[1:], chainIB[0]) +
    process.tracksSoA + process.beamSpotSoA + process.vertexSoA + process.unsortedOfflinePrimaryVerticesAlpaka + sum(chainAlpaka[1:], chainAlpaka[0]))


outputCommands = cms.untracked.vstring(list(process.FEVTDEBUGHLTEventContent.outputCommands) + [
    'keep *_*InBlocks*_*_REVTX', 'keep *_*PrimaryVerticesAlpaka*_*_REVTX', 'keep *_tofPIDAlpaka_*_REVTX',
    'drop *_offlinePrimaryVertices_*_RECO',
    'drop *_offlinePrimaryVertices4D_*_RECO',
    'drop *_offlinePrimaryVertices4DWithBS_*_RECO',
    'drop *_offlinePrimaryVerticesWithBS_*_RECO',
    'drop *_TriggerResults_*_RECO',
    'drop *_trackTimeValueMapProducer_generalTracksConfigurableFlatResolutionModel_RECO',
    'drop *_trackTimeValueMapProducer_generalTracksConfigurableFlatResolutionModelResolution_RECO',
    'drop *_trackTimeValueMapProducer_generalTracksPerfectResolutionModel_RECO',
    'drop *_trackTimeValueMapProducer_generalTracksPerfectResolutionModelResolution_RECO',
    'drop *_tofPID_probK_RECO',
    'drop *_tofPID_probP_RECO',
    'drop *_tofPID_probPi_RECO',
    'drop *_tofPID_sigmat0_RECO',
    'drop *_tofPID_sigmat0safe_RECO',
    'drop *_tofPID_t0_RECO',
    'drop *_tofPID_t0safe_RECO',
    'drop *_ak4CaloJetsForTrk_*_RECO',
    'drop *_inclusiveSecondaryVertices_*_RECO',
    'drop *_generalV0Candidates_Kshort_RECO',
    'drop *_generalV0Candidates_Lambda_RECO',
])
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(dataTier = cms.untracked.string('GEN-SIM-RECO'), filterName = cms.untracked.string('')),
    fileName = cms.untracked.string('file:revtx_step3_alpaka.root'),
    outputCommands = outputCommands,
    splitLevel = cms.untracked.int32(0),
)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

from Validation.Performance.TimeMemorySummary import customiseWithTimeMemorySummary
process = customiseWithTimeMemorySummary(process)
