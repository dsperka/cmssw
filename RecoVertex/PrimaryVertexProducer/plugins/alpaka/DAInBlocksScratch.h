#ifndef RecoVertex_PrimaryVertexProducer_plugins_alpaka_DAInBlocksScratch_h
#define RecoVertex_PrimaryVertexProducer_plugins_alpaka_DAInBlocksScratch_h

#include <Eigen/Core>
#include <Eigen/Dense>

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

// Scratch space of the DA-in-blocks clusterizer, per track and per vertex slot. It is only used by the clusterizer kernels
// (the track scratch has the same slots as tracksInBlocks, the vertex scratch the same as the vertex collection) and is
// neither put in the event nor copied to the host: the persistent TrackForVertexSoA and VertexSoA carry only what the
// fitter and the later modules read.
namespace reco {
  // Capacity of the per-track vertex vectors (maximum number of vertices, in total, that a track can be compared against)
  constexpr int maxTotalVertexForDeterministicAnnealing = 1024;

  using TrackToVertex = Eigen::Vector<float, maxTotalVertexForDeterministicAnnealing>;

  // Capacity of the per-vertex track lists (maximum number of tracks a single vertex can hold)
  constexpr int maxTracksPerVertexForDeterministicAnnealing = 1024;

  using VertexToTrack = Eigen::Vector<float, maxTracksPerVertexForDeterministicAnnealing>;
  using VertexToTrackInt = Eigen::Vector<int, maxTracksPerVertexForDeterministicAnnealing>;

  GENERATE_SOA_LAYOUT(TrackForVertexScratchSoALayout,
                      // columns: one value per track
                      SOA_COLUMN(float, sum_Z),  // partition function of the track
                      SOA_COLUMN(int, kmin),     // range of the vertices the track is compared with: kmin to kmax-1
                      SOA_COLUMN(int, kmax),
                      SOA_COLUMN(float, aux1),
                      SOA_COLUMN(float, aux2),

                      // Track-vertex association
                      SOA_EIGEN_COLUMN(TrackToVertex, vert_sw),
                      SOA_EIGEN_COLUMN(TrackToVertex, vert_se),
                      SOA_EIGEN_COLUMN(TrackToVertex, vert_swz),
                      SOA_EIGEN_COLUMN(TrackToVertex, vert_swE))

  using TrackForVertexScratchSoA = TrackForVertexScratchSoALayout<>;

  // Per-vertex-slot scratch, in the same slots as VertexSoA
  GENERATE_SOA_LAYOUT(VertexScratchSoALayout,
                      SOA_COLUMN(float, rho),  // regularisation (prior) weight of the vertex
                      SOA_COLUMN(float, aux1),
                      SOA_COLUMN(float, aux2),
                      SOA_COLUMN(float, sw),
                      SOA_COLUMN(float, swE),

                      // Dense track list of the vertex, filled by the clusterizer and read by the fitter. The compact
                      // version of it that is put in the event is VertexTrackSoA
                      SOA_EIGEN_COLUMN(VertexToTrackInt, track_id),
                      SOA_EIGEN_COLUMN(VertexToTrack, track_weight))

  using VertexScratchSoA = VertexScratchSoALayout<>;
}  // namespace reco

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using TrackForVertexScratchCollection = PortableCollection<::reco::TrackForVertexScratchSoA>;
  using VertexScratchCollection = PortableCollection<::reco::VertexScratchSoA>;
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoVertex_PrimaryVertexProducer_plugins_alpaka_DAInBlocksScratch_h
