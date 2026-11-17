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
  ALPAKA_FN_ACC __attribute__((noinline)) static void rejectOutliers(const TAcc& acc,
                                                                     TrackForVertexDeviceCollection::View tracks,
                                                                     TrackForVertexScratchCollection::View scratch,
                                                                     VertexDeviceCollection::View vertices,
                                                                     VertexScratchCollection::View vertexScratch,
                                                                     DAInBlocksClusterParameters const& cParams,
                                                                     double& osumtkwt,
                                                                     double& beta,
                                                                     int trackBlockSize) {
    // Treat outliers, either low quality vertex, or those with very far away tracks
    int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
    double rho0 = 0.0;
    if (cParams.dzCutOff > 0) {
      rho0 = vertices[blockIdx].nV() > 1 ? 1. / vertices[blockIdx].nV() : 1.;
      for (int rhoindex = 0; rhoindex < 5; rhoindex++) {
        update(acc,
               tracks,
               scratch,
               vertices,
               vertexScratch,
               cParams,
               osumtkwt,
               beta,
               rhoindex * rho0 / 5.,
               false,
               trackBlockSize);
      }
    }  // end if
    thermalize(
        acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, cParams.delta_lowT, rho0, trackBlockSize);
    int nprev = vertices[blockIdx].nV();
    merge(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
    while (nprev != vertices[blockIdx].nV()) {
      set_vtx_range(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
      update(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, rho0, false, trackBlockSize);
      nprev = vertices[blockIdx].nV();
      merge(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
    }
    while (beta < 1. / cParams.Tpurge) {
      alpaka::syncBlockThreads(acc);
      if (once_per_block(acc)) {
        beta = alpaka::math::min(acc, beta / cParams.coolingFactor, 1. / cParams.Tpurge);
      }
      alpaka::syncBlockThreads(acc);
      thermalize(acc,
                 tracks,
                 scratch,
                 vertices,
                 vertexScratch,
                 cParams,
                 osumtkwt,
                 beta,
                 cParams.delta_lowT,
                 rho0,
                 trackBlockSize);
    }
    alpaka::syncBlockThreads(acc);
    nprev = vertices[blockIdx].nV();
    purge(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, rho0, trackBlockSize);
    while (nprev != vertices[blockIdx].nV()) {
      thermalize(acc,
                 tracks,
                 scratch,
                 vertices,
                 vertexScratch,
                 cParams,
                 osumtkwt,
                 beta,
                 cParams.delta_lowT,
                 rho0,
                 trackBlockSize);
      nprev = vertices[blockIdx].nV();
      purge(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, rho0, trackBlockSize);
      alpaka::syncBlockThreads(acc);
    }
    while (beta < 1. / cParams.Tstop) {
      alpaka::syncBlockThreads(acc);
      if (once_per_block(acc)) {
        beta = alpaka::math::min(acc, beta / cParams.coolingFactor, 1. / cParams.Tstop);
      }
      alpaka::syncBlockThreads(acc);
      thermalize(acc,
                 tracks,
                 scratch,
                 vertices,
                 vertexScratch,
                 cParams,
                 osumtkwt,
                 beta,
                 cParams.delta_lowT,
                 rho0,
                 trackBlockSize);
    }
    alpaka::syncBlockThreads(acc);
    // The last track to vertex assignment!
    set_vtx_range(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
  }  // rejectOutliers

  class RejectOutliersKernel {
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
        printf("[DAInBlocksClusterizerAlgoRejectOutliers::operator()] Start outlier rejection for block %i\n",
               blockIdx);
        printf("[DAInBlocksClusterizerAlgoRejectOutliers::operator()] Parameter, trackBlockSize %i\n", trackBlockSize);
      }
#endif
      // After splitting we might get some candidates that are very low quality/have very far away tracks
      rejectOutliers(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgoRejectOutliers::operator()] End outlier rejection for block %i\n", blockIdx);
        printf("[STAGE-CMP] reject_outliers blockIdx=%d nV=%d\n", blockIdx, vertices[blockIdx].nV());
      }
#endif
      alpaka::syncBlockThreads(acc);
    }
  };  // class kernel

  void DAInBlocksClusterizerAlgo::reject_outliers(Queue& queue,
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
                        RejectOutliersKernel{},
                        deviceTrack.view(),
                        scratch.view(),
                        deviceVertex.view(),
                        vertexScratch.view(),
                        cParams,
                        beta_.data(),
                        osumtkwt_.data(),
                        blockSize);
  }  // DAInBlocksClusterizerAlgo::reject_outliers
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
