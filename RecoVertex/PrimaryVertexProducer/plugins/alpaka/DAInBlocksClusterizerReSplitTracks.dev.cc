#include "DAInBlocksClusterizerAlgo.h"
#include "DAInBlocksClusterizerPrimitives.h"
#include "DAInBlocksClusterizerMergeSplit.h"
#include "DAInBlocksClusterizerThermalize.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using namespace cms::alpakatools;
  using namespace daInBlocks;
  //////////////////////
  // Device functions //
  //////////////////////

  template <alpaka::concepts::Acc TAcc>
  ALPAKA_FN_ACC __attribute__((noinline)) static void reSplitTracks(const TAcc& acc,
                                                                    TrackForVertexDeviceCollection::View tracks,
                                                                    TrackForVertexScratchCollection::View scratch,
                                                                    VertexDeviceCollection::View vertices,
                                                                    VertexScratchCollection::View vertexScratch,
                                                                    DAInBlocksClusterParameters const& cParams,
                                                                    double& osumtkwt,
                                                                    double& beta,
                                                                    int trackBlockSize) {
    // Last splitting at the minimal temperature which is a bit more permissive
    int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
    int ntry = 0;
    double threshold = 1.0;
    int nprev = vertices[blockIdx].nV();
    split(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, threshold, trackBlockSize);
    while (nprev != vertices[blockIdx].nV() && (ntry++ < 10)) {
      thermalize(acc,
                 tracks,
                 scratch,
                 vertices,
                 vertexScratch,
                 cParams,
                 osumtkwt,
                 beta,
                 cParams.delta_highT,
                 0.0,
                 trackBlockSize);
      nprev = vertices[blockIdx].nV();
      merge(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
      while (nprev != vertices[blockIdx].nV()) {
        nprev = vertices[blockIdx].nV();
        update(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, 0.0, false, trackBlockSize);
        merge(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
      }
      threshold *= 1.1;
      split(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, threshold, trackBlockSize);
    }
  }

  class ReSplitTracksKernel {
  public:
    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC void operator()(const TAcc& acc,
                                  TrackForVertexDeviceCollection::View tracks,
                                  TrackForVertexScratchCollection::View scratch,
                                  VertexDeviceCollection::View vertices,
                                  VertexScratchCollection::View vertexScratch,
                                  DAInBlocksClusterParameters const cParams,
                                  double* beta_,
                                  double* osumtkwt_,
                                  int trackBlockSize) const {
      int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];  // Block number inside grid

      double& beta = alpaka::declareSharedVar<double, __COUNTER__>(acc);
      double& osumtkwt = alpaka::declareSharedVar<double, __COUNTER__>(acc);

      if (once_per_block(acc)) {
        beta = beta_[blockIdx];
        osumtkwt = osumtkwt_[blockIdx];
      }
      alpaka::syncBlockThreads(acc);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerReSplitTracks::operator()] BlockIdx %i, beta=%1.3f \n", blockIdx, beta);
      }
#endif
      // And split those with tension
      reSplitTracks(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerReSplitTracks::operator()] BlockIdx %i, end\n", blockIdx);
        printf("[STAGE-CMP] resplit_tracks blockIdx=%d nV=%d\n", blockIdx, vertices[blockIdx].nV());
      }
#endif
      alpaka::syncBlockThreads(acc);
    }
  };  // class kernel

  void DAInBlocksClusterizerAlgo::resplit_tracks(Queue& queue,
                                                 TrackForVertexDeviceCollection& deviceTrack,
                                                 TrackForVertexScratchCollection& scratch,
                                                 VertexDeviceCollection& deviceVertex,
                                                 VertexScratchCollection& vertexScratch,
                                                 DAInBlocksClusterParameters const& cParams,
                                                 int32_t nBlocks,
                                                 int32_t blockSize) {
    const int blocks = nBlocks;  // one block of blockSize threads per block of tracks
    alpaka::exec<Acc1D>(queue,
                        make_workdiv<Acc1D>(blocks, blockSize),
                        ReSplitTracksKernel{},
                        deviceTrack.view(),
                        scratch.view(),
                        deviceVertex.view(),
                        vertexScratch.view(),
                        cParams,
                        beta_.data(),
                        osumtkwt_.data(),
                        blockSize);
  }  // DAInBlocksClusterizerAlgo::resplit_tracks
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
