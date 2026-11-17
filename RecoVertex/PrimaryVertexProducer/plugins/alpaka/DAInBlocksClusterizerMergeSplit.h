#ifndef RecoVertex_PrimaryVertexProducer_plugins_alpaka_DAInBlocksClusterizerMergeSplit_h
#define RecoVertex_PrimaryVertexProducer_plugins_alpaka_DAInBlocksClusterizerMergeSplit_h

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
    ALPAKA_FN_ACC __attribute__((noinline)) inline void merge(const TAcc& acc,
                                                              TrackForVertexDeviceCollection::View tracks,
                                                              TrackForVertexScratchCollection::View scratch,
                                                              VertexDeviceCollection::View vertices,
                                                              VertexScratchCollection::View vertexScratch,
                                                              DAInBlocksClusterParameters const& cParams,
                                                              double& osumtkwt,
                                                              double& beta,
                                                              int trackBlockSize) {
      // If two vertex are too close together, merge them
      int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      int maxVerticesPerBlock = (int)1024 / alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      int nprev = vertices[blockIdx].nV();
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgo::merge()] BlockIdx %i, start merging \n", blockIdx);
      }
#endif
      if (nprev < 2)
        return;
      for (int ivertexB : independent_group_elements(acc, vertices[blockIdx].nV() - 1)) {
        int ivertexO = maxVerticesPerBlock * blockIdx + ivertexB;
        int ivertex = vertices[ivertexO].order();
        int ivertexnext = vertices[ivertexO + 1].order();
        vertexScratch[ivertex].aux1() = alpaka::math::abs(acc, vertices[ivertex].z() - vertices[ivertexnext].z());
      }
      alpaka::syncBlockThreads(acc);
      // Sorter things
      auto& critical_dist = alpaka::declareSharedVar<float[128], __COUNTER__>(acc);
      auto& critical_index = alpaka::declareSharedVar<int[128], __COUNTER__>(acc);
      int& ncritical = alpaka::declareSharedVar<int, __COUNTER__>(acc);
      if (once_per_block(acc)) {
        // Only thread 0 runs this, so step by 1: striding by blockSize here (as if all
        // threads shared the work) skipped most vertices on real GPU.
        ncritical = 0;
        for (int ivertexO = maxVerticesPerBlock * blockIdx;
             ivertexO < maxVerticesPerBlock * blockIdx + vertices[blockIdx].nV() - 1;
             ivertexO++) {
          int ivertex = vertices[ivertexO].order();
          if (vertexScratch[ivertex].aux1() < cParams.zmerge) {  // i.e., if we are to split the vertex
            critical_dist[ncritical] = alpaka::math::abs(acc, vertexScratch[ivertex].aux1());
            critical_index[ncritical] = ivertexO;
            ncritical++;
            if (ncritical > 128)
              break;
          }
        }
      }  // end once_per_block
      alpaka::syncBlockThreads(acc);
      if (ncritical == 0)
        return;
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgo::merge()] BlockIdx %i, %i vertices to be merged\n", blockIdx, ncritical);
      }
#endif
      int& ikO = alpaka::declareSharedVar<int, __COUNTER__>(acc);
      for (int sortO = 0; sortO < ncritical; ++sortO) {
        // All threads are running the same code, to know when to exit
        if (ncritical == 0 || maxVerticesPerBlock == nprev)
          return;
        // Only one thread picks the next vertex to merge: doing this redundantly on every
        // thread raced on critical_dist[ikO] (confirmed via compute-sanitizer racecheck).
        if (once_per_block(acc)) {
          ikO = 0;
          double minVal = 999999.;
          for (int sort1 = 0; sort1 < ncritical; ++sort1) {
            if (critical_dist[sort1] < minVal) {
              minVal = critical_dist[sort1];
              ikO = sort1;
            }
          }
          critical_dist[ikO] = 9999999.;
        }
        alpaka::syncBlockThreads(acc);
        int ivertexO = critical_index[ikO];
        int ivertex = vertices[ivertexO].order();
        int ivertexnext = blockIdx * maxVerticesPerBlock + nprev - 1;
        // A little bit of safety here. First is needed to avoid reading the -1 entry of vertices->order. Second is only needed if we go over 1024 vertices
        if (ivertexO < blockIdx * maxVerticesPerBlock + nprev - 1)
          ivertexnext = vertices[ivertexO + 1].order();
        alpaka::syncBlockThreads(acc);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
        if (once_per_block(acc)) {
          printf("[DAInBlocksClusterizerAlgo::merge()] BlockIdx %i, merge vertex %i into vertex %i\n",
                 blockIdx,
                 ivertex,
                 ivertexnext);
        }
#endif

        if (once_per_block(acc)) {
          vertices[ivertex].isGood() = false;  // Not deleting all the info, just disable it!
          double rho = vertexScratch[ivertex].rho() + vertexScratch[ivertexnext].rho();
          if (rho > 1.e-40) {
            vertices[ivertexnext].z() = (vertexScratch[ivertex].rho() * vertices[ivertex].z() +
                                         vertexScratch[ivertexnext].rho() * vertices[ivertexnext].z()) /
                                        rho;
          } else {
            vertices[ivertexnext].z() = 0.5 * (vertices[ivertex].z() + vertices[ivertexnext].z());
          }
          vertexScratch[ivertexnext].rho() = rho;
          vertexScratch[ivertexnext].sw() += vertexScratch[ivertex].sw();
          for (int ivertexOO = ivertexO; ivertexOO < maxVerticesPerBlock * blockIdx + nprev; ++ivertexOO) {
            vertices[ivertexOO].order() = vertices[ivertexOO + 1].order();
          }
          vertices[blockIdx].nV() = vertices[blockIdx].nV() - 1;
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
          if (once_per_block(acc)) {
            printf("[DAInBlocksClusterizerAlgo::merge()] BlockIdx %i, merged vertex %i with z=%1.3f,rho=%1.3f\n",
                   blockIdx,
                   ivertexnext,
                   vertices[ivertexnext].z(),
                   vertexScratch[ivertexnext].rho());
          }
#endif
        }  // end once_per_block
        alpaka::syncBlockThreads(acc);
        for (int resort = 0; resort < ncritical; ++resort) {
          if (critical_index[resort] > ivertexO)
            critical_index[resort]--;  // critical_index refers to the original vertices->order, so it needs to be updated
        }
        nprev = vertices[blockIdx].nV();  // And to the counter of previous vertices
        for (int itrackO : independent_group_elements(acc, trackBlockSize)) {
          int itrack = itrackO + blockIdx * trackBlockSize;
          if (scratch[itrack].kmax() > ivertexO + 1)
            scratch[itrack].kmax()--;
          if ((scratch[itrack].kmin() > ivertexO) || ((scratch[itrack].kmax() < (scratch[itrack].kmin() + 1)) &&
                                                      (scratch[itrack].kmin() > maxVerticesPerBlock * blockIdx)))
            scratch[itrack].kmin()--;
        }
        break;
      }
      alpaka::syncBlockThreads(acc);
      set_vtx_range(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
    }

    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC __attribute__((noinline)) inline void split(const TAcc& acc,
                                                              TrackForVertexDeviceCollection::View tracks,
                                                              TrackForVertexScratchCollection::View scratch,
                                                              VertexDeviceCollection::View vertices,
                                                              VertexScratchCollection::View vertexScratch,
                                                              DAInBlocksClusterParameters const& cParams,
                                                              double& osumtkwt,
                                                              double& beta,
                                                              double threshold,
                                                              int trackBlockSize) {
      int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      int maxVerticesPerBlock = (int)1024 / alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      update(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, 0.0, true, trackBlockSize);
      double epsilon = 1e-3;
      int nprev = vertices[blockIdx].nV();
      // Set critical T for all vertices
      for (int ivertexB : independent_group_elements(acc, vertices[blockIdx].nV())) {
        int ivertexO = maxVerticesPerBlock * blockIdx + ivertexB;
        int ivertex = vertices[ivertexO].order();
        double Tc = 2 * vertexScratch[ivertex].swE() / vertexScratch[ivertex].sw();
        vertexScratch[ivertex].aux1() = Tc;
      }
      alpaka::syncBlockThreads(acc);
      auto& critical_temp = alpaka::declareSharedVar<float[128], __COUNTER__>(acc);
      auto& critical_index = alpaka::declareSharedVar<int[128], __COUNTER__>(acc);
      int& ncritical = alpaka::declareSharedVar<int, __COUNTER__>(acc);
      // Information for the vertex splitting properties
      double& p1 = alpaka::declareSharedVar<double, __COUNTER__>(acc);
      double& p2 = alpaka::declareSharedVar<double, __COUNTER__>(acc);
      double& z1 = alpaka::declareSharedVar<double, __COUNTER__>(acc);
      double& z2 = alpaka::declareSharedVar<double, __COUNTER__>(acc);
      double& w1 = alpaka::declareSharedVar<double, __COUNTER__>(acc);
      double& w2 = alpaka::declareSharedVar<double, __COUNTER__>(acc);

      if (once_per_block(acc)) {
        // Only thread 0 runs this, so step by 1: striding by blockSize here (as if all
        // threads shared the work) skipped most vertices on real GPU.
        ncritical = 0;
        for (int ivertexO = maxVerticesPerBlock * blockIdx;
             ivertexO < maxVerticesPerBlock * blockIdx + vertices[blockIdx].nV();
             ivertexO++) {
          int ivertex = vertices[ivertexO].order();
          if (vertexScratch[ivertex].aux1() * beta > threshold) {
            critical_temp[ncritical] = alpaka::math::abs(acc, vertexScratch[ivertex].aux1());
            critical_index[ncritical] = ivertexO;
            ncritical++;
            if (ncritical > 128)
              break;
          }
        }
      }  // end once_per_block
      alpaka::syncBlockThreads(acc);
      if (ncritical == 0 || maxVerticesPerBlock == nprev)
        return;
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgo::split()] BlockIdx %i, split %i vertices\n", blockIdx, ncritical);
      }
#endif
      int& ikO = alpaka::declareSharedVar<int, __COUNTER__>(acc);
      for (int sortO = 0; sortO < ncritical; ++sortO) {  // All threads are running the same code, to know when to exit
        if (ncritical == 0 || maxVerticesPerBlock == nprev)
          return;
        // Only one thread picks the next vertex to split: doing this redundantly on every
        // thread raced on critical_temp[ikO] (confirmed via compute-sanitizer racecheck).
        if (once_per_block(acc)) {
          ikO = 0;
          double maxVal = -1.;
          for (int sort1 = 0; sort1 < ncritical; ++sort1) {
            if (critical_temp[sort1] > maxVal) {
              maxVal = critical_temp[sort1];
              ikO = sort1;
            }
          }
          critical_temp[ikO] = -1.;
        }
        alpaka::syncBlockThreads(acc);
        int ivertexO = critical_index[ikO];
        int ivertex = vertices[ivertexO].order();  // This will be splitted
        int ivertexprev = blockIdx * maxVerticesPerBlock;
        int ivertexnext = blockIdx * maxVerticesPerBlock + nprev - 1;
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
        if (once_per_block(acc)) {
          printf("[DAInBlocksClusterizerAlgo::split()] BlockIdx %i, splitting vertex %i\n", blockIdx, ivertex);
        }
#endif
        // Safety here. First is needed to avoid reading the -1 entry of vertices->order. Second in case we go over 511 vertices, but better keep it just in case
        if (ivertexO > blockIdx * maxVerticesPerBlock)
          ivertexprev = vertices[ivertexO - 1].order();
        if (ivertexO < blockIdx * maxVerticesPerBlock + nprev - 1)
          ivertexnext = vertices[ivertexO + 1].order();
        if (once_per_block(acc)) {
          p1 = 0.;
          p2 = 0.;
          z1 = 0.;
          z2 = 0.;
          w1 = 0.;
          w2 = 0.;
        }
        alpaka::syncBlockThreads(acc);
        for (int itrackO : independent_group_elements(acc, trackBlockSize)) {
          int itrack = itrackO + blockIdx * trackBlockSize;
          if (scratch[itrack].sum_Z() > 1e-40) {
            // winner-takes-all, usually overestimates splitting
            double tl = tracks[itrack].z() < vertices[ivertex].z() ? 1. : 0.;
            double tr = 1. - tl;
            // soften it, especially at low T
            double arg = (tracks[itrack].z() - vertices[ivertex].z()) *
                         alpaka::math::sqrt(acc, (beta)*tracks[itrack].oneoverdz2());
            if (alpaka::math::abs(acc, arg) < 20) {
              double t = alpaka::math::exp(acc, -arg);
              tl = t / (t + 1.);
              tr = 1 / (t + 1.);
            }
            // Recompute split vertex quantities
            double p =
                vertexScratch[ivertex].rho() * tracks[itrack].weight() *
                alpaka::math::exp(acc,
                                  -(beta) * (tracks[itrack].z() - vertices[ivertex].z()) *
                                      (tracks[itrack].z() - vertices[ivertex].z()) * tracks[itrack].oneoverdz2()) /
                scratch[itrack].sum_Z();
            double w = p * tracks[itrack].oneoverdz2();
            alpaka::atomicAdd(acc, &p1, p * tl, alpaka::hierarchy::Threads{});
            alpaka::atomicAdd(acc, &p2, p * tr, alpaka::hierarchy::Threads{});
            alpaka::atomicAdd(acc, &z1, w * tl * tracks[itrack].z(), alpaka::hierarchy::Threads{});
            alpaka::atomicAdd(acc, &z2, w * tr * tracks[itrack].z(), alpaka::hierarchy::Threads{});
            alpaka::atomicAdd(acc, &w1, w * tl, alpaka::hierarchy::Threads{});
            alpaka::atomicAdd(acc, &w2, w * tr, alpaka::hierarchy::Threads{});
          }
        }
        alpaka::syncBlockThreads(acc);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
        if (once_per_block(acc)) {
          printf(
              "[DAInBlocksClusterizerAlgo::split()] BlockIdx %i, for vertex %i, p1=%1.3f, p2=%1.3f, w1=%1.3f, "
              "w2=%1.3f, "
              "z1=%1.3f, "
              "z2=%1.3f\n",
              blockIdx,
              ivertex,
              p1,
              p2,
              w1,
              w2,
              z1,
              z2);
        }
#endif

        if (once_per_block(acc)) {
          // If one vertex is taking all the things, then set the others slightly off to help splitting
          if (w1 > 0) {
            z1 = z1 / w1;
          } else {
            z1 = vertices[ivertex].z() - epsilon;
          }
          if (w2 > 0) {
            z2 = z2 / w2;
          } else {
            z2 = vertices[ivertex].z() + epsilon;
          }
          // If there is not enough room, reduce split size
          if ((ivertexO > maxVerticesPerBlock * blockIdx) &&
              (z1 <
               (0.6 * vertices[ivertex].z() +
                0.4 *
                    vertices[ivertexprev]
                        .z()))) {  // First in the if is ivertexO, as we care on whether the vertex is the leftmost or rightmost
            z1 = 0.6 * vertices[ivertex].z() + 0.4 * vertices[ivertexprev].z();
          }
          if ((ivertexO < maxVerticesPerBlock * blockIdx + nprev - 1) &&
              (z2 > (0.6 * vertices[ivertex].z() + 0.4 * vertices[ivertexnext].z()))) {
            z2 = 0.6 * vertices[ivertex].z() + 0.4 * vertices[ivertexnext].z();
          }
        }  // end once_per_block
        alpaka::syncBlockThreads(acc);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
        if (once_per_block(acc)) {
          printf("[DAInBlocksClusterizerAlgo::split()] BlockIdx %i, vertex %i will split into z1=%1.3f, z2=%1.3f\n",
                 blockIdx,
                 ivertex,
                 z1,
                 z2);
        }
#endif
        // Split only if the two sub-clusters are significantly separated and their weights are finite (as in the C++
        // implementation): otherwise leave this vertex alone and go on with the next critical one
        bool& doSplit = alpaka::declareSharedVar<bool, __COUNTER__>(acc);
        if (once_per_block(acc)) {
          const double pk1 = p1 * vertexScratch[ivertex].rho() / (p1 + p2);
          const double pk2 = p2 * vertexScratch[ivertex].rho() / (p1 + p2);
          doSplit = ((z2 - z1) > epsilon) && alpaka::math::isfinite(acc, pk1) && alpaka::math::isfinite(acc, pk2);
        }
        alpaka::syncBlockThreads(acc);
        if (!doSplit)
          continue;
        const int noFreeSlot = 999999;  // marks 'no empty vertex slot found'
        int nnew = noFreeSlot;
        // Find the first empty index to save the vertex
        for (int icheck = maxVerticesPerBlock * blockIdx; icheck < maxVerticesPerBlock * (blockIdx + 1); icheck++) {
          if (not(vertices[icheck].isGood())) {
            nnew = icheck;
            break;
          }
        }
        if (nnew == noFreeSlot)
          break;  // No free slot left: leave the splitting loop.
                  // Every thread runs the free-slot search above itself (instead of one thread publishing the result)
        // and they all read the same vertex state, so all of them take this branch together and exit consistently.
        if (once_per_block(acc)) {
          double pk1 = p1 * vertexScratch[ivertex].rho() / (p1 + p2);
          double pk2 = p2 * vertexScratch[ivertex].rho() / (p1 + p2);
          vertices[ivertex].z() = z2;
          vertexScratch[ivertex].rho() = pk2;
          // Insert it into the first available slot
          vertices[nnew].z() = z1;
          vertexScratch[nnew].rho() = pk1;
          vertices[nnew].isGood() = true;
          // This is likely not needed as far as it is reset anytime we call update but better be safe in case we reenable a previously disable vertex (i.e. split in the memory where a merged one was)
          vertexScratch[nnew].sw() = 0.;
          vertexScratch[nnew].swE() = 0.;
          for (int ivnew = maxVerticesPerBlock * blockIdx + nprev; ivnew > ivertexO;
               ivnew--) {  // As we add a vertex, we update from the back downwards
            vertices[ivnew].order() = vertices[ivnew - 1].order();
          }
          vertices[ivertexO].order() = nnew;
          vertices[blockIdx].nV() += 1;
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
          if (once_per_block(acc)) {
            printf("[DAInBlocksClusterizerAlgo::split()] BlockIdx %i, vertex %i did split into indexes %i and %i\n",
                   blockIdx,
                   ivertex,
                   ivertex,
                   nnew);
            dump(acc, beta, vertices, vertexScratch);
          }
#endif
        }
        alpaka::syncBlockThreads(acc);
        // Now, update kmin/kmax for all tracks
        for (int itrackO : independent_group_elements(acc, trackBlockSize)) {
          int itrack = itrackO + blockIdx * trackBlockSize;
          if (scratch[itrack].kmin() > ivertexO)
            scratch[itrack].kmin()++;
          if ((scratch[itrack].kmax() >= ivertexO) || (scratch[itrack].kmax() == scratch[itrack].kmin()))
            scratch[itrack].kmax()++;
        }
        nprev = vertices[blockIdx].nV();
        if (once_per_block(acc)) {
          // If we did a splitting or old sorted list of vertex index is scrambled, so we need to fix it
          for (int resort = 0; resort < ncritical; ++resort) {
            if (critical_index[resort] > ivertexO)
              critical_index[resort]++;
          }
        }
        alpaka::syncBlockThreads(acc);
      }
      alpaka::syncBlockThreads(acc);
    }

    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC __attribute__((noinline)) inline void purge(const TAcc& acc,
                                                              TrackForVertexDeviceCollection::View tracks,
                                                              TrackForVertexScratchCollection::View scratch,
                                                              VertexDeviceCollection::View vertices,
                                                              VertexScratchCollection::View vertexScratch,
                                                              DAInBlocksClusterParameters const& cParams,
                                                              double& osumtkwt,
                                                              double& beta,
                                                              double rho0,
                                                              int trackBlockSize) {
      // Remove repetitive or low quality entries
      int blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      int maxVerticesPerBlock = (int)1024 / alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      if (vertices[blockIdx].nV() < 2)
        return;
      double eps = 1e-40;
      int nunique_min = 2;
      double rhoconst = rho0 * alpaka::math::exp(acc, -beta * (cParams.dzCutOff * cParams.dzCutOff));
      int nprev = vertices[blockIdx].nV();
      // Reassign
      set_vtx_range(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
      for (int ivertexB : independent_group_elements(acc, vertices[blockIdx].nV())) {
        int ivertexO = maxVerticesPerBlock * blockIdx + ivertexB;
        int ivertex = vertices[ivertexO].order();
        vertexScratch[ivertex].aux1() = 0;  // sum of track-vertex probabilities
        vertexScratch[ivertex].aux2() = 0;  // number of uniquely assigned tracks
      }
      alpaka::syncBlockThreads(acc);
      // Get quality of vertex in terms of #Tracks and sum of track probabilities
      for (int itrackO : independent_group_elements(acc, trackBlockSize)) {
        int itrack = itrackO + blockIdx * trackBlockSize;
        double track_aux1 = ((scratch[itrack].sum_Z() > eps) && (tracks[itrack].weight() > cParams.uniquetrkminp))
                                ? 1. / scratch[itrack].sum_Z()
                                : 0.;
        for (int ivertexO = scratch[itrack].kmin(); ivertexO < scratch[itrack].kmax(); ++ivertexO) {
          int ivertex = vertices[ivertexO].order();
          double ppcut =
              cParams.uniquetrkweight * vertexScratch[ivertex].rho() / (vertexScratch[ivertex].rho() + rhoconst);
          double track_vertex_aux1 = alpaka::math::exp(
              acc,
              -(beta)*tracks[itrack].oneoverdz2() *
                  ((tracks[itrack].z() - vertices[ivertex].z()) * (tracks[itrack].z() - vertices[ivertex].z())));
          float p = vertexScratch[ivertex].rho() * track_vertex_aux1 *
                    track_aux1;  // The whole track-vertex P_ij = rho_j*p_ij*p_i
          alpaka::atomicAdd(acc, &vertexScratch[ivertex].aux1(), p, alpaka::hierarchy::Threads{});
          if (p > ppcut) {
            alpaka::atomicAdd(acc, &vertexScratch[ivertex].aux2(), 1.f, alpaka::hierarchy::Threads{});
          }
        }
      }
      alpaka::syncBlockThreads(acc);
      // Find worst vertex to purge
      int& k0 = alpaka::declareSharedVar<int, __COUNTER__>(acc);

      if (once_per_block(acc)) {
        // Only thread 0 runs this, so step by 1 (same fix as in split()/merge()).
        double sumpmin = tracks.nT();
        k0 = maxVerticesPerBlock * blockIdx + nprev;
        for (int ivertexO = maxVerticesPerBlock * blockIdx;
             ivertexO < maxVerticesPerBlock * blockIdx + (int)vertices[blockIdx].nV();
             ivertexO++) {
          int ivertex = vertices[ivertexO].order();
          if ((vertexScratch[ivertex].aux2() < nunique_min) && (vertexScratch[ivertex].aux1() < sumpmin)) {
            sumpmin = vertexScratch[ivertex].aux1();
            k0 = ivertexO;
          }
        }  // end vertex for
        if (k0 != (int)(maxVerticesPerBlock * blockIdx + nprev)) {
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
          if (once_per_block(acc)) {
            printf(
                "[DAInBlocksClusterizerAlgo::purge()] BlockIdx %i, some vertices need purging. Will start purging \n",
                blockIdx);
          }
#endif
          int ivertexRemoved = vertices[k0].order();
          for (int ivertexOO = k0; ivertexOO < maxVerticesPerBlock * blockIdx + (int)nprev - 1; ++ivertexOO) {
            vertices[ivertexOO].order() =
                vertices[ivertexOO + 1].order();  // Update vertex order taking out the purged one
          }
          vertices[blockIdx].nV()--;
          vertices[ivertexRemoved].isGood() = false;
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_CLUSTERIZERALGO
          if (once_per_block(acc)) {
            printf("[DAInBlocksClusterizerAlgo::purge()] BlockIdx %i, vertex %i purged\n", blockIdx, k0);
          }
#endif
        }
      }  // end once_per_block
      alpaka::syncBlockThreads(acc);
      if (k0 != (int)(maxVerticesPerBlock * blockIdx + (int)nprev)) {
        for (int itrackO : independent_group_elements(acc, trackBlockSize)) {
          int itrack = itrackO + blockIdx * trackBlockSize;
          if (scratch[itrack].kmax() > k0)
            scratch[itrack].kmax()--;
          if ((scratch[itrack].kmin() > k0) || ((scratch[itrack].kmax() < (scratch[itrack].kmin() + 1)) &&
                                                (scratch[itrack].kmin() > (int)(maxVerticesPerBlock * blockIdx))))
            scratch[itrack].kmin()--;
        }
      }  // end if
      alpaka::syncBlockThreads(acc);
      if (nprev != vertices[blockIdx].nV()) {
        set_vtx_range(acc, tracks, scratch, vertices, vertexScratch, cParams, osumtkwt, beta, trackBlockSize);
      }
    }

  }  // namespace daInBlocks

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoVertex_PrimaryVertexProducer_plugins_alpaka_DAInBlocksClusterizerMergeSplit_h
