#ifndef RecoVertex_PrimaryVertexProducer_plugins_alpaka_TrackForVertexScratch_h
#define RecoVertex_PrimaryVertexProducer_plugins_alpaka_TrackForVertexScratch_h

#include <Eigen/Core>
#include <Eigen/Dense>

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

// Per-track scratch space of the DA-in-blocks clusterizer. It is only used by the clusterizer kernels, in the same order of
// tracks as tracksInBlocks, and is neither put in the event nor copied to the host: the persistent TrackForVertexSoA
// carries only the track quantities and the track-vertex assignment that later modules read.
namespace reco {
  // Capacity of the per-track vertex vectors (maximum number of vertices, in total, that a track can be compared against)
  constexpr int maxTotalVertexForDeterministicAnnealing = 1024;

  using TrackToVertex = Eigen::Vector<float, maxTotalVertexForDeterministicAnnealing>;

  GENERATE_SOA_LAYOUT(TrackForVertexScratchSoALayout,
                      // columns: one value per track
                      SOA_COLUMN(float, aux1),
                      SOA_COLUMN(float, aux2),

                      // Track-vertex association
                      SOA_EIGEN_COLUMN(TrackToVertex, vert_sw),
                      SOA_EIGEN_COLUMN(TrackToVertex, vert_se),
                      SOA_EIGEN_COLUMN(TrackToVertex, vert_swz),
                      SOA_EIGEN_COLUMN(TrackToVertex, vert_swE))

  using TrackForVertexScratchSoA = TrackForVertexScratchSoALayout<>;
}  // namespace reco

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using TrackForVertexScratchCollection = PortableCollection<::reco::TrackForVertexScratchSoA>;
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoVertex_PrimaryVertexProducer_plugins_alpaka_TrackForVertexScratch_h
