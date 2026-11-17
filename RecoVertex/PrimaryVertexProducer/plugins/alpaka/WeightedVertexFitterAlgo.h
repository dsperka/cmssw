#ifndef RecoVertex_PrimaryVertexProducer_plugins_alpaka_WeightedVertexFitterAlgo_h
#define RecoVertex_PrimaryVertexProducer_plugins_alpaka_WeightedVertexFitterAlgo_h

#include "DAInBlocksScratch.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexTrackDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/TrackForVertexDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // the beam spot quantities the fit uses, passed to the kernel by value (no device copy)
  struct BeamSpotForFit {
    float x, y;
    float beamWidthX, beamWidthY;
  };

  struct FitterParameters {
    // In principle a single one, but in case we end up wanting to add more configuration to the fiter
    bool useBeamSpotConstraint;
  };

  class WeightedVertexFitterAlgo {
  public:
    WeightedVertexFitterAlgo(Queue& queue, FitterParameters fPar);
    void fit(Queue& queue,
             const TrackForVertexDeviceCollection& deviceTrack,
             VertexDeviceCollection& deviceVertex,
             VertexScratchCollection& vertexScratch,
             VertexTrackDeviceCollection& vertexTracks,
             BeamSpotForFit const& beamSpot);

  private:
    bool useBeamSpotConstraint_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoVertex_PrimaryVertexProducer_plugins_alpaka_WeightedVertexFitterAlgo_h
