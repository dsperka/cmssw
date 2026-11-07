#ifndef RecoVertex_PrimaryVertexProducer_plugins_alpaka_DAInBlocksClusterizerThermalize_h
#define RecoVertex_PrimaryVertexProducer_plugins_alpaka_DAInBlocksClusterizerThermalize_h

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "DataFormats/VertexSoA/interface/alpaka/VertexDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/TrackForVertexDeviceCollection.h"

#include "DAInBlocksClusterizerAlgo.h"
#include "DAInBlocksClusterizerPrimitives.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using namespace cms::alpakatools;

  namespace daInBlocks {

    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC __attribute__((noinline)) inline void thermalize(const TAcc& acc,
                                                                   TrackForVertexDeviceCollection::View tracks,
                                                                   VertexDeviceCollection::View vertices,
                                                                   DAInBlocksClusterParameters const& cParams,
                                                                   double& osumtkwt,
                                                                   double& beta,
                                                                   double delta_highT,
                                                                   double rho0,
                                                                   int trackBlockSize) {
      // At a fixed temperature, iterate vertex position update until stable
      int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      int maxVerticesPerBlock = (int)1024 / alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      int niter = 0;
      float zrange_min = 0.01;  // Hard coded as it is in CPU
      float delta_max = cParams.delta_lowT;
      if (cParams.convergence_mode == 0) {
        delta_max = delta_highT;
      } else if (cParams.convergence_mode == 1) {
        delta_max = cParams.delta_lowT / alpaka::math::sqrt(acc, alpaka::math::max(acc, beta, 1.0));
      }
      int maxIterations = 1000;
      set_vtx_range(acc, tracks, vertices, cParams, osumtkwt, beta, trackBlockSize);
      // Accumulator of variations
      float delta_sum_range = 0;
      while (niter++ < maxIterations) {
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
        if (once_per_block(acc)) {
          printf("[DAInBlocksClusterizerAlgo::thermalize()] BlockIdx %i, thermalize at beta=%1.3f, iteration %i\n",
                 blockIdx,
                 beta,
                 niter);
        }
#endif
        update(acc, tracks, vertices, cParams, osumtkwt, beta, rho0, false, trackBlockSize);
        float dmax = 0.;
        for (int ivertexO = maxVerticesPerBlock * blockIdx;
             ivertexO < maxVerticesPerBlock * blockIdx + vertices[blockIdx].nV();
             ivertexO++) {
          int ivertex = vertices[ivertexO].order();
          if (alpaka::math::abs(acc, vertices[ivertex].aux1()) > dmax)
            dmax = alpaka::math::abs(acc, vertices[ivertex].aux1());
        }
        delta_sum_range += dmax;
        // If a vertex moved too much we reassign
        if (delta_sum_range > zrange_min && dmax > zrange_min) {
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
          if (once_per_block(acc)) {
            printf(
                "[DAInBlocksClusterizerAlgo::thermalize()] BlockIdx %i, thermalize at beta=%1.3f, iteration %i. Found "
                "delta_sum_range=%1.3f, dmax=%1.3f, will redo track-vertex assignament\n",
                blockIdx,
                beta,
                niter,
                delta_sum_range,
                dmax);
          }
#endif
          set_vtx_range(acc, tracks, vertices, cParams, osumtkwt, beta, trackBlockSize);
          delta_sum_range = 0.;
        }
        if (dmax < delta_max) {
          update(acc, tracks, vertices, cParams, osumtkwt, beta, 0.0, true, trackBlockSize);

#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
          if (once_per_block(acc)) {
            printf(
                "[DAInBlocksClusterizerAlgo::thermalize()] BlockIdx %i, thermalize at beta=%1.3f, iteration %i. Found "
                "delta_sum_range=%1.3f, dmax=%1.3f, all vertices stable enough to stop thermalizing\n",
                blockIdx,
                beta,
                niter,
                delta_sum_range,
                dmax);
            dump(acc, beta, vertices);
          }
#endif
          break;
        }
      }  // end while
    }  // thermalize

  }  // namespace daInBlocks

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoVertex_PrimaryVertexProducer_plugins_alpaka_DAInBlocksClusterizerThermalize_h
