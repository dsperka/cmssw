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
                 minNdof = cms.double(0.0), useBeamConstraint = cms.bool(False), maxDistanceToBeam = cms.double(1.0),
                 vertexTimeParameters = cms.PSet(algorithm = cms.string('fromTracksPID'))),
        cms.PSet(label = cms.string("WithBS"), algorithm = cms.string("WeightedMeanFitter"), chi2cutoff = cms.double(2.5),
                 minNdof = cms.double(0.0), useBeamConstraint = cms.bool(True), maxDistanceToBeam = cms.double(1.0),
                 vertexTimeParameters = cms.PSet(algorithm = cms.string('fromTracksPID')))))

# Alpaka: portable tracks -> DA in blocks on the device -> reco::Vertex
process.tracksSoA = cms.EDProducer("PortableTrackSoAProducer@alpaka",
    TrackLabel = cms.InputTag("generalTracks"),
    BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
    TkFilterParameters = cms.PSet(
        maxNormalizedChi2 = cms.double(10.0), minPixelLayersWithHits = cms.int32(2), minSiliconLayersWithHits = cms.int32(5),
        maxD0Significance = cms.double(4.0), maxD0Error = cms.double(1.0), maxDzError = cms.double(1.0), minPt = cms.double(0.0),
        maxEta = cms.double(4.0), trackQuality = cms.string("any"), vertexSize = cms.double(0.006), d0CutOff = cms.double(3.)))
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
    soaVertex = cms.InputTag("vertexSoA"), srcTrack = cms.InputTag("generalTracks"),
    # vertex time as for the standard 3D collection (fromTracksPID)
    vertexTimeParameters = cms.PSet(algorithm = cms.string("fromTracksPID")))

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


# Slim event content: only what the MTD vertex validation reads (about 66 MB/event instead of 170 MB/event with the full
# FEVTDEBUGHLT content, and the job runs about 1.7 times faster with 6 threads). The validation of a file written with this
# list gives the same histograms as the one of the full content (checked on 10 events, 3191 histograms, differences < 1e-15).
# Not droppable: TrackingParticles/TrackingVertexs (truth, 38 MB/event), GenParticles (referenced by the truth), the Gsf
# electron products (mtdEleIsoValid resolves their references).
outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_generalTracks_*_*', 'keep *_trackExtenderWithMTD_*_*', 'keep *_mtdTrackQualityMVA_*_*',
    'keep *_mtdClusters_*_*', 'keep *_mtdRecHits_*_*', 'keep *_mtdUncalibratedRecHits_*_*', 'keep *_mtdTrackingRecHits_*_*',
    'keep *_tofPID*_*_*', 'keep *_*rimaryVertices*_*_*', 'keep *_offlineBeamSpot_*_*', 'keep *_generatorSmeared_*_*',
    'keep *Gsf*_*_*_*', 'keep *_electronGsfTracks*_*_*', 'keep *_gedGsfElectrons_*_*', 'keep *_gedGsfElectronCores_*_*',
    'keep *_ecalDrivenGsfElectronsHGC_*_*', 'keep *_ecalDrivenGsfElectronCoresHGC_*_*', 'keep *_ecalDrivenGsfElectronCores_*_*',
    'keep *_particleFlowSuperClusterECAL_*_*', 'keep *_particleFlowSuperClusterHGCal_*_*', 'keep *_particleFlowSuperClusterHGCalFromTICL_*_*',
    'keep *_particleFlowClusterECAL_*_*', 'keep *_particleFlowClusterHGCal_*_*', 'keep *_particleFlowClusterHGCalFromTICL_*_*',
    'keep *_hgcalMergeLayerClusters_*_*',
    'keep recoGenParticles_*_*_*', 'keep *_genParticles_*_*', 'keep edmHepMCProduct_*_*_*', 'keep GenEventInfoProduct_*_*_*',
    'keep *_mix_MergedTrackTruth_*', 'keep *_mix_FTLBarrel_*', 'keep *_mix_FTLEndcap_*', 'keep *_mix_MergedMtdTruth*_*',
    'keep *_g4SimHits_FastTimerHits*_*',
    'keep *_trackingParticleRecoTrackAsssociation_*_*', 'keep *_mtdRecoClusterToSimLayerClusterAssociation_*_*',
    'keep *_mtdSimLayerClusterToTPAssociation_*_*',
    'keep *_*_*_REVTX',
    # debug SoA collections of the Alpaka and GNN chains, not read by the validation
    'drop *_vertexSoA_*_*', 'drop *_gnnVertexProducer_*_*', 'drop *PortableHostCollection*_*_*_*',
)
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(dataTier = cms.untracked.string('GEN-SIM-RECO'), filterName = cms.untracked.string('')),
    fileName = cms.untracked.string('file:revtx_step3_alpaka.root'),
    outputCommands = outputCommands,
    splitLevel = cms.untracked.int32(0),
)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

from Validation.Performance.TimeMemorySummary import customiseWithTimeMemorySummary
process = customiseWithTimeMemorySummary(process)
