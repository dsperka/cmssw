#ifndef RecoVertex_PrimaryVertexProducer_plugins_alpaka_DAInBlocksClusterizerPrimitives_h
#define RecoVertex_PrimaryVertexProducer_plugins_alpaka_DAInBlocksClusterizerPrimitives_h

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "DataFormats/VertexSoA/interface/alpaka/VertexDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/TrackForVertexDeviceCollection.h"

#include "DAInBlocksClusterizerAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using namespace cms::alpakatools;

  namespace daInBlocks {

#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC inline void dump(const TAcc& acc,
                                   double& beta,
                                   VertexDeviceCollection::View vertices,
                                   VertexScratchCollection::View vertexScratch) {
      int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      int maxVerticesPerBlock =
          (int)1024 / alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(
                          acc)[0u];  // Max vertices size is 1024 over number of blocks in the kernel
      printf("[DAInBlocksClusterizerAlgo::dump()] Block Idx %i with nV %i at beta %1.8f \n",
             blockIdx,
             vertices[blockIdx].nV(),
             beta);
      for (int ivertexO = maxVerticesPerBlock * blockIdx;
           ivertexO < maxVerticesPerBlock * blockIdx + vertices[blockIdx].nV();
           ivertexO += 1) {
        int ivertex = vertices[ivertexO].order();
        printf(
            "[DAInBlocksClusterizerAlgo::dump()] -- Block Idx %i, vertex %i in order %i: "
            "z=%1.5f,swE=%1.5f,sw=%1.5f,pk=%1.5f\n",
            blockIdx,
            ivertex,
            ivertexO,
            vertices[ivertex].z(),
            vertexScratch[ivertex].swE(),
            vertexScratch[ivertex].sw(),
            vertexScratch[ivertex].rho());
      }
    }
#endif  // DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO

    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC __attribute__((noinline)) inline void set_vtx_range(const TAcc& acc,
                                                                      TrackForVertexDeviceCollection::View tracks,
                                                                      TrackForVertexScratchCollection::View scratch,
                                                                      VertexDeviceCollection::View vertices,
                                                                      VertexScratchCollection::View vertexScratch,
                                                                      DAInBlocksClusterParameters const& cParams,
                                                                      double& osumtkwt,
                                                                      double& beta,
                                                                      int trackBlockSize) {
      // These updates the range of vertices associated to each track through the kmin/kmax variables
      int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];  // Block number inside grid
      int maxVerticesPerBlock = (int)1024 / alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(
                                                acc)[0u];  // Max vertices size is 1024 over number of blocks in grid
      float zrange_min = 0.1;                              // Hard coded as in CPU version
      for (int itrackO : independent_group_elements(acc, trackBlockSize)) {
        int itrack = itrackO + blockIdx * trackBlockSize;
        // Based on current temperature (regularization term) and track position uncertainty, only keep relevant vertices
        float zrange = alpaka::math::max(
            acc, cParams.zrange / alpaka::math::sqrt(acc, (beta)*tracks[itrack].oneoverdz2()), (double)zrange_min);
        float zmin = tracks[itrack].z() - zrange;
        // First the lower bound
        int kmin = alpaka::math::min(
            acc, (int)(maxVerticesPerBlock * blockIdx) + vertices.nV(blockIdx) - 1, scratch[itrack].kmin());
        // If the vertex position in z is bigger than the minimum, go down through all vertices position until finding one that is too far
        if (vertices[vertices[kmin].order()].z() > zmin) {
          // i.e., while we find another vertex within range that is before the previous initial step
          while ((kmin > maxVerticesPerBlock * blockIdx) && ((vertices[vertices[kmin - 1].order()].z()) > zmin)) {
            kmin--;
          }
        } else {
          while ((kmin < (maxVerticesPerBlock * blockIdx + (int)(vertices[blockIdx].nV()) - 1)) &&
                 ((vertices[vertices[kmin].order()].z()) < zmin)) {
            kmin++;
          }
        }

        // And now do the same for the upper bound
        float zmax = tracks[itrack].z() + zrange;
        int kmax =
            alpaka::math::max(acc,
                              0,
                              alpaka::math::min(acc,
                                                maxVerticesPerBlock * blockIdx + (int)(vertices[blockIdx].nV()) - 1,
                                                (int)(scratch[itrack].kmax()) - 1));
        // For corner cases in which we purged the first vertex, thus not properly updating kmax during purging
        while (vertices[kmax].order() == emptyVertexOrder)
          kmax++;
        if (vertices[vertices[kmax].order()].z() < zmax) {
          // As long as we have more vertex above kmax but within z range, we can add them to the collection, keep going
          while ((kmax < (maxVerticesPerBlock * blockIdx + (int)(vertices[blockIdx].nV()) - 1)) &&
                 ((vertices[vertices[kmax + 1].order()].z()) < zmax)) {
            kmax++;
          }
        } else {
          while ((kmax > maxVerticesPerBlock * blockIdx) && (vertices[vertices[kmax].order()].z() > zmax)) {
            kmax--;
          }
        }
        if (kmin <= kmax) {
          scratch[itrack].kmin() = (int)kmin;
          scratch[itrack].kmax() = (int)kmax + 1;
        } else {  // Track goes to the most extreme cases if no associated one is found
          scratch[itrack].kmin() =
              (int)alpaka::math::max(acc, maxVerticesPerBlock * blockIdx, (int)alpaka::math::min(acc, kmin, kmax));
          scratch[itrack].kmax() =
              (int)alpaka::math::min(acc,
                                     (maxVerticesPerBlock * blockIdx) + (int)vertices[blockIdx].nV(),
                                     (int)alpaka::math::max(acc, kmin, kmax) + 1);
        }
      }  //end for
      alpaka::syncBlockThreads(acc);
    }

    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC __attribute__((noinline)) inline void update(const TAcc& acc,
                                                               TrackForVertexDeviceCollection::View tracks,
                                                               TrackForVertexScratchCollection::View scratch,
                                                               VertexDeviceCollection::View vertices,
                                                               VertexScratchCollection::View vertexScratch,
                                                               DAInBlocksClusterParameters const& cParams,
                                                               double& osumtkwt,
                                                               double& beta,
                                                               double rho0,
                                                               bool updateTc,
                                                               int trackBlockSize) {
      // Main function that updates the annealing parameters on each T step, computes all partition functions and so on
      int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];  // Block number inside grid
      int maxVerticesPerBlock = (int)1024 / alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(
                                                acc)[0u];  // Max vertices size is 1024 over number of blocks in grid
      // Initial partition function, really only used on the outlier rejection step to penalize
      double Zinit = rho0 * alpaka::math::exp(acc, -(beta)*cParams.dzCutOff * cParams.dzCutOff);
      for (int itrackO : independent_group_elements(acc, trackBlockSize)) {
        int itrack = itrackO + blockIdx * trackBlockSize;
        double botrack_dz2 = -(beta)*tracks[itrack].oneoverdz2();
        scratch[itrack].sum_Z() = Zinit;
        for (int ivertexO = scratch[itrack].kmin(); ivertexO < scratch[itrack].kmax(); ++ivertexO) {
          int ivertex = vertices[ivertexO].order();
          if (not(vertices[ivertex].isGood()))
            continue;
          // Z_t = sum_v pho_v * e^{-beta*(z_t-z_v)/dz^2}, partition function of the track
          scratch[itrack].sum_Z() += vertexScratch[ivertex].rho() *
                                     alpaka::math::exp(acc,
                                                       botrack_dz2 * (tracks[itrack].z() - vertices[ivertex].z()) *
                                                           (tracks[itrack].z() - vertices[ivertex].z()));
        }  //end vertex for
        if (not(alpaka::math::isfinite(acc, scratch[itrack].sum_Z())))
          scratch[itrack].sum_Z() = 0;
        if (scratch[itrack].sum_Z() > 1e-40) {  // If non-zero then the track has a non-trivial assignment to a vertex
          double sumw = tracks[itrack].weight() / (double)scratch[itrack].sum_Z();
          for (int ivertexO = scratch[itrack].kmin(); ivertexO < scratch[itrack].kmax(); ++ivertexO) {
            int ivertex = vertices[ivertexO].order();
            if (not(vertices[ivertex].isGood()))
              continue;
            double vert_exparg = botrack_dz2 * (tracks[itrack].z() - vertices[ivertex].z()) *
                                 (tracks[itrack].z() - vertices[ivertex].z());  // -beta*(z_t-z_v)/dz^2
            double vert_exp = alpaka::math::exp(acc, vert_exparg);              // e^{-beta*(z_t-z_v)/dz^2}
            scratch[itrack].vert_se()[ivertex] =
                vert_exp * sumw;  // From partition of track to contribution of track to vertex partition
            double w = vertexScratch[ivertex].rho() * vert_exp * sumw * tracks[itrack].oneoverdz2();
            scratch[itrack].vert_sw()[ivertex] = w;
            scratch[itrack].vert_swz()[ivertex] = w * tracks[itrack].z();
            if (updateTc) {
              scratch[itrack].vert_swE()[ivertex] = -w * vert_exparg / (beta);
            } else {
              scratch[itrack].vert_swE()[ivertex] = 0;
            }
          }  //end vertex for
        }  //end if
      }  //end track for
      alpaka::syncBlockThreads(acc);
      // After the track-vertex matrix assignment, we need to add up across vertices. This time, we use one thread per vertex
      for (int ivertexB : independent_group_elements(acc, vertices[blockIdx].nV())) {
        int ivertexO = maxVerticesPerBlock * blockIdx + ivertexB;
        int ivertex = vertices[ivertexO].order();
        float se = 0.;
        float sw = 0.;
        float swz = 0.;
        float swE = 0.;
        for (int itrack = blockIdx * trackBlockSize; itrack < (blockIdx + 1) * trackBlockSize; itrack++) {
          if ((ivertexO >= scratch[itrack].kmin()) && (ivertexO < scratch[itrack].kmax())) {
            se += scratch[itrack].vert_se()[ivertex];
            sw += scratch[itrack].vert_sw()[ivertex];
            swz += scratch[itrack].vert_swz()[ivertex];
            swE += scratch[itrack].vert_swE()[ivertex];
          }
        }
        if (sw > 0) {
          double znew = swz / sw;
          vertexScratch[ivertex].aux1() = znew - vertices[ivertex].z();
          vertices[ivertex].z() = znew;
        } else {
          vertexScratch[ivertex].aux1() = 0.;
        }
        vertexScratch[ivertex].rho() = vertexScratch[ivertex].rho() * se * osumtkwt;
        vertexScratch[ivertex].sw() = sw;
        if (updateTc)
          vertexScratch[ivertex].swE() = swE;
      }  // end vertex for
      alpaka::syncBlockThreads(acc);
    }  //end update

  }  // namespace daInBlocks

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoVertex_PrimaryVertexProducer_plugins_alpaka_DAInBlocksClusterizerPrimitives_h
