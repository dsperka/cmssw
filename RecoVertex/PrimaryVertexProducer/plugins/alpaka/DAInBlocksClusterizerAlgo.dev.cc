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
  ALPAKA_FN_ACC static void initialize(const TAcc& acc,
                                       TrackForVertexDeviceCollection::View tracks,
                                       TrackForVertexScratchCollection::View scratch,
                                       VertexDeviceCollection::View vertices,
                                       VertexScratchCollection::View vertexScratch,
                                       DAInBlocksClusterParameters const& cParams,
                                       int trackBlockSize) {
    // Initialize all vertices as empty, a single vertex in each block will be initialized with all tracks associated to it
    int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
    int maxVerticesPerBlock = (int)1024 / alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u];
    vertices[blockIdx].nV() = 1;
    for (int ivertexI : independent_group_elements(
             acc,
             maxVerticesPerBlock)) {  // Initialize vertices in parallel in the block. Note that a block of threads should always be restricted to oeprations within maxVerticesPerBlock*blockIdx and maxVerticesPerBlock*(blockIdx+1)-1 to avoid running conditions
      int ivertex = ivertexI + maxVerticesPerBlock * blockIdx;
      vertexScratch[ivertex].sw() = 0.;
      vertexScratch[ivertex].swE() = 0.;
      vertices[ivertex].z() = 0.;
      vertexScratch[ivertex].rho() = 0.;
      vertices[ivertex].isGood() = false;
      vertices[ivertex].order() = emptyVertexOrder;
      if (ivertex ==
          maxVerticesPerBlock *
              blockIdx) {  // Set up the initial single vetex containing everything which should only happen for the first vertex
        vertexScratch[ivertex].rho() = 1.;
        vertices[ivertex].order() = maxVerticesPerBlock * blockIdx;
        vertices[ivertex].isGood() = true;
      }
    }  // end for
    alpaka::syncBlockThreads(acc);
    for (int itrackO : independent_group_elements(acc, trackBlockSize)) {
      int itrack = itrackO + blockIdx * trackBlockSize;
      // Tracks are associated to vertex in list kmin, kmin+1,... kmax-1, so this just assign all tracks to the vertex we just created
      scratch[itrack].kmin() = maxVerticesPerBlock * blockIdx;
      scratch[itrack].kmax() = maxVerticesPerBlock * blockIdx + 1;
    }
    alpaka::syncBlockThreads(acc);
  }

  template <alpaka::concepts::Acc TAcc>
  ALPAKA_FN_ACC static void getBeta0(const TAcc& acc,
                                     TrackForVertexDeviceCollection::View tracks,
                                     TrackForVertexScratchCollection::View scratch,
                                     VertexDeviceCollection::View vertices,
                                     VertexScratchCollection::View vertexScratch,
                                     DAInBlocksClusterParameters const& cParams,
                                     double& beta,
                                     int trackBlockSize) {
    // Computes first critical temperature
    int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
    int maxVerticesPerBlock = (int)1024 / alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u];
    for (int itrackO : independent_group_elements(acc, trackBlockSize)) {
      int itrack = itrackO + blockIdx * trackBlockSize;
      if (not(tracks[itrack].isGood()))
        continue;
      scratch[itrack].aux1() = tracks[itrack].weight() * tracks[itrack].oneoverdz2();
      scratch[itrack].aux2() = tracks[itrack].weight() * tracks[itrack].oneoverdz2() * tracks[itrack].z();
    }
    // Initial vertex position
    alpaka::syncBlockThreads(acc);
    float& wnew = alpaka::declareSharedVar<float, __COUNTER__>(acc);
    float& znew = alpaka::declareSharedVar<float, __COUNTER__>(acc);
    if (once_per_block(acc)) {
      wnew = 0.;
      znew = 0.;
    }
    alpaka::syncBlockThreads(acc);
    for (int itrackO : independent_group_elements(acc, trackBlockSize)) {
      int itrack = itrackO + blockIdx * trackBlockSize;
      if (not(tracks[itrack].isGood()))
        continue;
      alpaka::atomicAdd(acc, &wnew, scratch[itrack].aux1(), alpaka::hierarchy::Threads{});
      alpaka::atomicAdd(acc, &znew, scratch[itrack].aux2(), alpaka::hierarchy::Threads{});
    }
    alpaka::syncBlockThreads(acc);
    if (once_per_block(acc)) {
      vertices[maxVerticesPerBlock * blockIdx].z() = znew / wnew;
      znew = 0.;
    }
    alpaka::syncBlockThreads(acc);
    // Now do a chi-2 like of all tracks and save it again in znew
    for (int itrackO : independent_group_elements(acc, trackBlockSize)) {
      int itrack = itrackO + blockIdx * trackBlockSize;
      if (not(tracks[itrack].isGood()))
        continue;
      scratch[itrack].aux2() =
          scratch[itrack].aux1() * (vertices[maxVerticesPerBlock * blockIdx].z() - tracks[itrack].z()) *
          (vertices[maxVerticesPerBlock * blockIdx].z() - tracks[itrack].z()) * tracks[itrack].oneoverdz2();
      alpaka::atomicAdd(acc, &znew, scratch[itrack].aux2(), alpaka::hierarchy::Threads{});
    }
    alpaka::syncBlockThreads(acc);
    if (once_per_block(acc)) {
      beta = 2 * znew / wnew;  // 1/beta_C, or T_C
      if (beta > cParams.Tmin) {
        int coolingsteps =
            1 - int(alpaka::math::log(acc, beta / cParams.Tmin) /
                    alpaka::math::log(
                        acc, cParams.coolingFactor));  // A tricky conversion to round the number of cooling steps
        beta = alpaka::math::pow(acc, cParams.coolingFactor, coolingsteps) / cParams.Tmin;
      } else {
        beta = cParams.coolingFactor / cParams.Tmin;
      }
    }
    alpaka::syncBlockThreads(acc);
  }

  template <alpaka::concepts::Acc TAcc>
  ALPAKA_FN_ACC static void coolingWhileSplitting(const TAcc& acc,
                                                  TrackForVertexDeviceCollection::View tracks,
                                                  TrackForVertexScratchCollection::View scratch,
                                                  VertexDeviceCollection::View vertices,
                                                  VertexScratchCollection::View vertexScratch,
                                                  DAInBlocksClusterParameters const& cParams,
                                                  double& osumtkwt,
                                                  double& beta,
                                                  int trackBlockSize) {
    // Perform cooling of the deterministic annealing
    int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
    double betafreeze = (1. / cParams.Tmin) * alpaka::math::sqrt(acc, cParams.coolingFactor);
    while (beta < betafreeze) {
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf(
            "[DAInBlocksClusterizerAlgo::coolingWhileSplitting()] BlockIdx %i, current beta=%1.8f\n", blockIdx, beta);
      }
#endif
      int nprev = vertices[blockIdx].nV();
      merge(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
      while (nprev != vertices[blockIdx].nV()) {
        // If we are here, we merged before, keep merging until stable
        nprev = vertices[blockIdx].nV();
        update(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, 0.0, false, trackBlockSize);
        merge(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
      }  // end while after merging
      split(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, 1.0, trackBlockSize);
      if (once_per_block(acc)) {
        // Cool down
        beta = beta / cParams.coolingFactor;
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
                 cParams.delta_highT,
                 0.0,
                 trackBlockSize);
    }
    set_vtx_range(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
    update(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, 0.0, false, trackBlockSize);
  }  // end coolingWhileSplitting

  template <alpaka::concepts::Acc TAcc>
  ALPAKA_FN_ACC static void reMergeTracks(const TAcc& acc,
                                          TrackForVertexDeviceCollection::View tracks,
                                          TrackForVertexScratchCollection::View scratch,
                                          VertexDeviceCollection::View vertices,
                                          VertexScratchCollection::View vertexScratch,
                                          DAInBlocksClusterParameters const& cParams,
                                          double& osumtkwt,
                                          double& beta,
                                          int trackBlockSize) {
    // After the cooling, we merge any closeby vertices
    int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
    int nprev = vertices[blockIdx].nV();
    merge(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
    while (nprev != vertices[blockIdx].nV()) {
      // Keep merging until stable
      set_vtx_range(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
      update(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, 0.0, false, trackBlockSize);
      nprev = vertices[blockIdx].nV();
      merge(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
    }  // end while
  }  // end reMergeTracks

  class DAInBlocksClusterizerKernel {
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
      // Core of the clusterization algorithm
      // Produces set of clusters for input set of block-overlapped tracks
      // tracks contains input track parameters and includes the track-vertex assignment modified during this kernel
      // vertices is filled up by this kernel with protocluster properties
      // beta_ and osumtkwt_ are used to pass the final values of beta and osumtkwt on each block to the next kernel
      int blockSize = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u];
      int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgo::operator()] Start clustering block %i\n", blockIdx);
        printf("[DAInBlocksClusterizerAlgo::operator()] Parameters blockSize %i, trackBlockSize %i\n",
               blockSize,
               trackBlockSize);
      }
#endif
      double& beta =
          alpaka::declareSharedVar<double, __COUNTER__>(acc);  // 1/T in the annealing loop, shared in the block
      double& osumtkwt = alpaka::declareSharedVar<double, __COUNTER__>(
          acc);  // Sum of all track weights for normalization of probabilities, shared in the block
      if (once_per_block(acc)) {
        osumtkwt = 0;
      }
      alpaka::syncBlockThreads(acc);
      for (int itrackO : independent_group_elements(acc, trackBlockSize)) {
        int itrack = itrackO + blockIdx * trackBlockSize;
        if (not(tracks[itrack].isGood()))
          continue;
        double temp_weight = static_cast<double>(tracks[itrack].weight());
        alpaka::atomicAdd(acc, &osumtkwt, temp_weight, alpaka::hierarchy::Threads{});
      }
      alpaka::syncBlockThreads(acc);
      if (once_per_block(acc)) {
        osumtkwt = 1. / osumtkwt;
      }
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgo::operator()] BlockIdx %i, osumtkwt=%1.3f \n", blockIdx, osumtkwt);
      }
#endif
      alpaka::syncBlockThreads(acc);
      // In each block, initialize to a single vertex with all tracks
      initialize(acc, tracks, scratch, vertices, vertexScratch, cParams, trackBlockSize);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgo::operator()] BlockIdx %i, vertices initialized\n", blockIdx);
        printf("[STAGE-CMP] after_initialize blockIdx=%d nV=%d\n", blockIdx, vertices[blockIdx].nV());
      }
#endif
      alpaka::syncBlockThreads(acc);
      // First estimation of critical temperature
      getBeta0(acc, tracks, scratch, vertices, vertexScratch, cParams, beta, trackBlockSize);
      alpaka::syncBlockThreads(acc);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgo::operator()] BlockIdx %i, first estimation of TC beta=%1.8f \n",
               blockIdx,
               beta);
        printf("[STAGE-CMP] after_getBeta0 blockIdx=%d beta=%.10f nV=%d\n", blockIdx, beta, vertices[blockIdx].nV());
      }
#endif
      // Cool down to betamax with rho = 0.0 (no regularization term)
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
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgo::operator()] BlockIdx %i, first thermalization ended\n", blockIdx);
        printf("[STAGE-CMP] after_thermalize blockIdx=%d beta=%.10f nV=%d\n", blockIdx, beta, vertices[blockIdx].nV());
      }
#endif
      // Now the cooling loop
      coolingWhileSplitting(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgo::operator()] BlockIdx %i, cooling ended, T at stop beta=%1.8f\n",
               blockIdx,
               beta);
        printf("[STAGE-CMP] after_coolingWhileSplitting blockIdx=%d beta=%.10f nV=%d\n",
               blockIdx,
               beta,
               vertices[blockIdx].nV());
      }
#endif
      // After cooling, merge closeby vertices
      reMergeTracks(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgo::operator()] BlockIdx %i, merge, last merging step done\n", blockIdx);
        printf("[STAGE-CMP] after_reMergeTracks blockIdx=%d nV=%d\n", blockIdx, vertices[blockIdx].nV());
      }
#endif
      if (once_per_block(acc)) {
        beta_[blockIdx] = beta;
        osumtkwt_[blockIdx] = osumtkwt;
      }
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgo::operator()] BlockIdx %i end\n", blockIdx);
        printf("[STAGE-CMP] clusterize blockIdx=%d nV=%d\n", blockIdx, vertices[blockIdx].nV());
      }
#endif
    }
  };  // DAInBlocksClusterizerKernel

  DAInBlocksClusterizerAlgo::DAInBlocksClusterizerAlgo(Queue& queue, int32_t bSize)
      : beta_(cms::alpakatools::make_device_buffer<double[]>(queue, bSize)),
        osumtkwt_(cms::alpakatools::make_device_buffer<double[]>(queue, bSize)) {
    alpaka::memset(queue, beta_, 0.);
    alpaka::memset(queue, osumtkwt_, 0.);
  }

  void DAInBlocksClusterizerAlgo::clusterize(Queue& queue,
                                             TrackForVertexDeviceCollection& deviceTrack,
                                             TrackForVertexScratchCollection& scratch,
                                             VertexDeviceCollection& deviceVertex,
                                             VertexScratchCollection& vertexScratch,
                                             DAInBlocksClusterParameters const& cParams,
                                             int32_t nBlocks,
                                             int32_t blockSize) {
    alpaka::exec<Acc1D>(queue,
                        make_workdiv<Acc1D>(nBlocks, 64),
                        DAInBlocksClusterizerKernel{},
                        deviceTrack.view(),
                        scratch.view(),
                        deviceVertex.view(),
                        vertexScratch.view(),
                        cParams,
                        beta_.data(),
                        osumtkwt_.data(),
                        blockSize);
  }  // DAInBlocksClusterizerAlgo::clusterize
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
