import FWCore.ParameterSet.Config as cms
from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import offlinePrimaryVertices as _offlinePrimaryVertices

# Primary vertices with the Alpaka implementation of the deterministic annealing in blocks
# (see Configuration.ProcessModifiers.vertexInBlocksAlpaka_cff).
#
#   pvTracksSoA --\
#   pvBeamSpotSoA -+-> pvVertexSoA       -> unsortedOfflinePrimaryVerticesFromSoA        (collection without beam constraint)
#                  +-> pvVertexSoAWithBS -> unsortedOfflinePrimaryVerticesFromSoAWithBS   (collection with beam constraint)
#
# The two converted collections replace the "" and "WithBS" collections of unsortedOfflinePrimaryVertices.

pvTracksSoA = cms.EDProducer("PortableTrackSoAProducer@alpaka",
    TrackLabel = cms.InputTag("generalTracks"),
    BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
    # same track selection as the standard producer; vertexSize and d0CutOff are used to compute the track weights
    TkFilterParameters = _offlinePrimaryVertices.TkFilterParameters.clone(
        vertexSize = cms.double(0.006),
        d0CutOff = cms.double(3.)
    )
)
# higher eta cut for the phase 2 tracker, as in the standard producer
from Configuration.Eras.Modifier_phase2_tracker_cff import phase2_tracker
phase2_tracker.toModify(pvTracksSoA, TkFilterParameters = dict(maxEta = 4.0))

pvBeamSpotSoA = cms.EDProducer("BeamSpotDeviceProducer@alpaka",
    src = cms.InputTag("offlineBeamSpot")
)

pvVertexSoA = cms.EDProducer("PrimaryVertexProducerPortable@alpaka",
    TrackLabel = cms.InputTag("pvTracksSoA"),
    BeamSpotLabel = cms.InputTag("pvBeamSpotSoA"),
    blockOverlap = cms.double(0.50),
    blockSize = cms.int32(512),
    TkFitterParameters = cms.PSet(
        useBeamSpotConstraint = cms.bool(False)
    ),
    TkClusParameters = cms.PSet(
        coolingFactor = cms.double(0.6),   # moderate annealing speed
        zrange = cms.double(4.),           # consider only clusters within 4 sigma*sqrt(T) of a track
        delta_highT = cms.double(1.e-2),   # convergence requirement at high T
        delta_lowT = cms.double(1.e-3),    # convergence requirement at low T
        convergence_mode = cms.int32(0),   # 0 = two steps, 1 = dynamic with sqrt(T)
        Tmin = cms.double(2.0),            # end of vertex splitting
        Tpurge = cms.double(2.0),          # cleaning
        Tstop = cms.double(0.5),           # end of annealing
        vertexSize = cms.double(0.006),    # added in quadrature to track-z resolutions
        d0CutOff = cms.double(3.),         # downweight high IP tracks
        dzCutOff = cms.double(3.),         # outlier rejection after freeze-out (T<Tmin)
        zmerge = cms.double(1e-2),         # merge intermediate clusters separated by less than zmerge
        uniquetrkweight = cms.double(0.8), # require at least two tracks with this weight at T=Tpurge
        uniquetrkminp = cms.double(0.0),   # minimal a priori track weight for counting unique tracks
        mintrkweight = cms.double(0.25)    # minimal probability for a track to be assigned to a vertex
    )
)
pvVertexSoAWithBS = pvVertexSoA.clone(
    TkFitterParameters = dict(useBeamSpotConstraint = True)
)

# vertex time as for the standard producer: only with the timing layer
unsortedOfflinePrimaryVerticesFromSoA = cms.EDProducer("RecoVertexFromSoAProducer",
    soaVertex = cms.InputTag("pvVertexSoA"),
    srcTrack = cms.InputTag("generalTracks"),
    vertexTimeParameters = cms.PSet(algorithm = cms.string(""))
)
unsortedOfflinePrimaryVerticesFromSoAWithBS = unsortedOfflinePrimaryVerticesFromSoA.clone(
    soaVertex = "pvVertexSoAWithBS"
)
from Configuration.Eras.Modifier_phase2_timing_layer_cff import phase2_timing_layer
phase2_timing_layer.toModify(unsortedOfflinePrimaryVerticesFromSoA,
    vertexTimeParameters = dict(algorithm = 'fromTracksPID'))
phase2_timing_layer.toModify(unsortedOfflinePrimaryVerticesFromSoAWithBS,
    vertexTimeParameters = dict(algorithm = 'fromTracksPID'))
