import FWCore.ParameterSet.Config as cms

# This modifier replaces the primary vertex producer by the Alpaka implementation
# of the deterministic annealing in blocks (runs on the device selected by the
# accelerator configuration, CPU or GPU):
# - reco::Tracks -> portable track SoA, beam spot -> portable beam spot
# - clustering in blocks and weighted vertex fit on the device
# - conversion of the portable vertices back to reco::Vertex (with the vertex time, if the timing layer is present)
# Used together with the vertexInBlocks modifier name space of the C++ implementation, but independent of it.
vertexInBlocksAlpaka = cms.Modifier()
