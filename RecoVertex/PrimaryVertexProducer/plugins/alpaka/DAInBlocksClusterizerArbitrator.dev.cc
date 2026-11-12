#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaInterface/interface/radixSort.h"

#include "DAInBlocksClusterizerAlgo.h"

//#define DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ARBITRATOR 1

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using namespace cms::alpakatools;
  //////////////////////
  // Device functions //
  //////////////////////

  template <alpaka::concepts::Acc TAcc>
  ALPAKA_FN_ACC static void resortVerticesAndAssign(const TAcc& acc,
                                                    TrackForVertexDeviceCollection::View tracks,
                                                    TrackForVertexScratchCollection::View scratch,
                                                    VertexDeviceCollection::View vertices,
                                                    DAInBlocksClusterParameters const& cParams,
                                                    int32_t nBlocks) {
    // Multiblock vertex arbitration
    double beta = 1. / cParams.Tstop;
    const unsigned int maxVerticesInSoA = 1024;
    auto& z = alpaka::declareSharedVar<float[maxVerticesInSoA], __COUNTER__>(acc);
    auto& rho = alpaka::declareSharedVar<float[maxVerticesInSoA], __COUNTER__>(acc);
    alpaka::syncBlockThreads(acc);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ARBITRATOR
    if (once_per_block(acc)) {
      printf("[DAInBlocksClusterizerAlgoArbitrator::resortVerticesAndAssign()] Start reassignment\n");
    }
#endif

    for (unsigned int iv : uniform_elements(acc, maxVerticesInSoA)) {
      if (vertices[iv].isGood()) {
        if ((vertices[iv].rho() > 10000) || (alpaka::math::abs(acc, vertices[iv].z()) > 30)) {
          vertices[iv].isGood() = false;
        }
      }
    }

#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ARBITRATOR
    if (once_per_block(acc)) {
      printf("[DAInBlocksClusterizerAlgoArbitrator::resortVerticesAndAssign()] Reorder vertices\n");
    }
#endif

    alpaka::syncBlockThreads(acc);
    if (once_per_block(acc)) {
      // We only keep the ones that are good out of all of the ones in the collection
      unsigned int nTrueVertex = 0;
      // Each block b owns the slots [b*maxVerticesPerBlock, (b+1)*maxVerticesPerBlock) and its order() list is only
      // valid for the first nV(b) entries: the entries behind them are stale and point at vertices that were merged
      // or purged (or already listed), so they must not be read
      const int maxVerticesPerBlock = (int)maxVerticesInSoA / nBlocks;
      for (int b = 0; b < nBlocks; b++) {
        for (int i = 0; i < vertices[b].nV(); i++) {
          int iv = vertices[b * maxVerticesPerBlock + i].order();
          if ((iv < 0) || (iv >= (int)maxVerticesInSoA))
            continue;
          if (vertices[iv].isGood()) {
            z[nTrueVertex] = vertices[iv].z();
            rho[nTrueVertex] = vertices[iv].rho();
            nTrueVertex++;
          }
        }
      }
      vertices[0].nV() = nTrueVertex;
    }
    alpaka::syncBlockThreads(acc);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ARBITRATOR
    if (once_per_block(acc)) {
      printf(
          "[DAInBlocksClusterizerAlgoArbitrator::resortVerticesAndAssign()] Before vertices are reordered, nV %i "
          "left\n",
          vertices[0].nV());
    }
#endif

    auto& orderedIndices = alpaka::declareSharedVar<uint16_t[maxVerticesInSoA], __COUNTER__>(acc);
    auto& sws = alpaka::declareSharedVar<uint16_t[maxVerticesInSoA], __COUNTER__>(acc);

    cms::alpakatools::radixSort<Acc1D, float, 2>(acc, z, orderedIndices, sws, vertices[0].nV());
    alpaka::syncBlockThreads(acc);
    // copy sorted vertices back to the SoA. We restrict our usage to the first vertices[0].nV() entries of the SoA
    for (int ivtx : uniform_elements(acc, vertices[0].nV())) {
      vertices[ivtx].z() = z[ivtx];
      vertices[ivtx].rho() = rho[ivtx];
      vertices[ivtx].order() = orderedIndices[ivtx];
      vertices[ivtx].isGood() = true;
    }
    // And invalidate the remaining part we won't use anymore
    for (unsigned int ivtxO : uniform_elements(acc, maxVerticesInSoA - vertices[0].nV())) {
      unsigned int ivtx = ivtxO + vertices[0].nV();
      vertices[ivtx].isGood() = false;
    }
    alpaka::syncBlockThreads(acc);
    float zrange_min = 0.1;

#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ARBITRATOR
    if (once_per_block(acc)) {
      printf("[DAInBlocksClusterizerAlgoArbitrator::resortVerticesAndAssign()] Vertices are reordered, nV %i left\n",
             vertices[0].nV());
    }
#endif

    for (int itrack : uniform_elements(acc, tracks.nT())) {
      if (not(tracks[itrack].isGood()))
        continue;
      float zrange = alpaka::math::max(
          acc, cParams.zrange / alpaka::math::sqrt(acc, (beta)*tracks[itrack].oneoverdz2()), (double)zrange_min);
      float zmin = tracks[itrack].z() - zrange;
      int kmin = vertices[0].nV() - 1;
      if (kmin < 0)
        continue;
      if (vertices[vertices[kmin].order()].z() > zmin) {
        while ((kmin > 0) && (vertices[vertices[kmin - 1].order()].z() > zmin)) {
          // i.e., while we find another vertex within range that is before the previous initial step we select it as minimum
          kmin--;
        }
      } else {
        while ((kmin < vertices[0].nV()) && (vertices[vertices[kmin].order()].z() < zmin)) {
          kmin++;
        }
      }
      // Now the same for the upper bound
      float zmax = tracks[itrack].z() + zrange;
      int kmax = 0;
      if (vertices[vertices[kmax].order()].z() < zmax) {
        while ((kmax < vertices[0].nV() - 1) && (vertices[vertices[kmax + 1].order()].z() < zmax)) {
          // As long as we have more vertex above kmax but within z range, we can add them to the collection, keep going
          kmax++;
        }
      } else {
        while ((kmax > 0) && (vertices[vertices[kmax].order()].z() > zmax)) {
          kmax--;
        }
      }
      if (kmin <= kmax) {
        scratch[itrack].kmin() = kmin;
        scratch[itrack].kmax() = kmax + 1;
      } else {
        scratch[itrack].kmin() = alpaka::math::max(acc, 0, alpaka::math::min(acc, kmin, kmax));
        scratch[itrack].kmax() = alpaka::math::min(acc, vertices[0].nV(), alpaka::math::max(acc, kmin, kmax) + 1);
      }
    }
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ARBITRATOR
    if (once_per_block(acc)) {
      printf(
          "[DAInBlocksClusterizerAlgoArbitrator::resortVerticesAndAssign()] Track-vertex range assignment finished \n");
    }
#endif

    alpaka::syncBlockThreads(acc);
    double rho0 = vertices[0].nV() > 1 ? 1. / vertices[0].nV() : 1.;
    double z_sum_init = rho0 * alpaka::math::exp(acc, -(beta)*cParams.dzCutOff * cParams.dzCutOff);
    for (int itrack : uniform_elements(acc, tracks.nT())) {
      int kmin = scratch[itrack].kmin();
      int kmax = scratch[itrack].kmax();
      double p_max = -1;
      const int noVertexAssigned = 10000;  // track has no vertex above mintrkweight
      int iMax = noVertexAssigned;
      double sum_Z = z_sum_init;
      for (auto k = kmin; k < kmax; k++) {
        double v_exp = alpaka::math::exp(
            acc,
            -(beta)*alpaka::math::pow(acc, tracks[itrack].z() - vertices[vertices[k].order()].z(), 2) *
                tracks[itrack].oneoverdz2());
        sum_Z += vertices[vertices[k].order()].rho() * v_exp;
      }
      double invZ = sum_Z > 1e-40 ? 1. / sum_Z : 0.0;
      // Univocally assign a track to the highest probability vertex
      for (auto k = kmin; k < kmax; k++) {
        float v_exp = alpaka::math::exp(
            acc,
            -(beta)*alpaka::math::pow(acc, tracks[itrack].z() - vertices[vertices[k].order()].z(), 2) *
                tracks[itrack].oneoverdz2());
        float p = vertices[vertices[k].order()].rho() * v_exp * invZ;
        if (p > p_max && p > cParams.mintrkweight) {
          // assign  track i -> vertex k (hard assignment unles we configure a non standard mintrweight_
          p_max = p;
          iMax = k;
        }
      }
      // Finally register that itrack is assigned to vertex iMax
      scratch[itrack].kmin() = iMax;
      scratch[itrack].kmax() = iMax + 1;
    }
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ARBITRATOR
    if (once_per_block(acc)) {
      printf("[DAInBlocksClusterizerAlgoArbitrator::resortVerticesAndAssign()] Updated vertex masses\n");
    }
#endif

    alpaka::syncBlockThreads(acc);
  }  //resortVerticesAndAssign

  template <alpaka::concepts::Acc TAcc>
  ALPAKA_FN_ACC static void finalizeVertices(const TAcc& acc,
                                             TrackForVertexDeviceCollection::View tracks,
                                             TrackForVertexScratchCollection::View scratch,
                                             VertexDeviceCollection::View vertices,
                                             DAInBlocksClusterParameters const& cParams) {
    // Per-candidate z uncertainty (inverse-variance combination of its
    // assigned tracks' oneoverdz2, i.e. 1/sigma_z^2), computed in parallel
    // below alongside each vertex's track list, then read serially by the
    // merge decision further down -- that decision has an inherent
    // left-to-right dependency so it can't itself be parallelized, but
    // computing each candidate's sigma doesn't depend on that chain at all.
    constexpr unsigned int maxVerticesInSoA = 1024;
    auto& sigma_z = alpaka::declareSharedVar<double[maxVerticesInSoA], __COUNTER__>(acc);
    // First put the tracks in vertex SoA
    for (int k : uniform_elements(acc, vertices[0].nV())) {
      int ivertex = vertices[k].order();
      vertices[ivertex].ntracks() = 0;
      for (int itrack = 0; itrack < tracks.nT(); itrack += 1) {
        if (not(tracks[itrack].isGood())) {
          continue;
        }
        int ivtxFromTk = scratch[itrack].kmin();
        if (ivtxFromTk == k) {
          bool isNew = true;
          for (int ivtrack = 0; ivtrack < vertices[ivertex].ntracks(); ivtrack++) {
            if (tracks[itrack].ttIndex() == tracks[vertices[ivertex].track_id()[ivtrack]].ttIndex())
              isNew = false;
          }
          if (!isNew) {
            continue;
          }
          // Here we are saving the itrack (i.e. index in the device collection) instead tracks[itrack].ttIndex(); (index of the input reco track collection) so we can use the same device collection in the fitter
          vertices[ivertex].track_id()[vertices[ivertex].ntracks()] = itrack;
          vertices[ivertex].track_weight()[vertices[ivertex].ntracks()] = 1.;
          vertices[ivertex].ntracks()++;
        }
      }
      double sum_invvar = 0.;
      for (int it = 0; it < vertices[ivertex].ntracks(); it++) {
        sum_invvar += tracks[vertices[ivertex].track_id()[it]].oneoverdz2();
      }
      sigma_z[ivertex] = (sum_invvar > 0.) ? alpaka::math::sqrt(acc, 1. / sum_invvar) : 0.;
      // Candidates with fewer than 2 tracks are not dropped here: as in the C++ implementation they take part in the
      // merge chain below and only the merged cluster has to hold more than one track
      vertices[ivertex].x() = 0;
      vertices[ivertex].y() = 0;
    }
    alpaka::syncBlockThreads(acc);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ARBITRATOR
    if (once_per_block(acc)) {
      printf("[DAInBlocksClusterizerAlgoArbitrator::finalizeVertices()] Before cleanup %i vertices\n",
             vertices[0].nV());
    }
#endif
    if (once_per_block(acc)) {
      // Merge left to right exactly as DAClusterizerInZ_vect::vertices_in_blocks does: consecutive candidates k and
      // k+1 (in z order, empty ones included) belong to the same cluster if they are closer than the larger of
      // 2*vertexSize and 3 combined sigma; the tracks of a cluster are collected in its first vertex and the cluster
      // is kept if it holds more than one track
      const int nCand = vertices[0].nV();
      int head = -1;
      for (int k = 0; k < nCand; k++) {
        const int thisVertex = (int)vertices[k].order();
        if (head < 0) {
          head = thisVertex;
        } else {
          for (int itrack = 0; itrack < vertices[thisVertex].ntracks(); itrack++) {
            vertices[head].track_id()[vertices[head].ntracks()] = vertices[thisVertex].track_id()[itrack];
            vertices[head].track_weight()[vertices[head].ntracks()] = 1.;
            vertices[head].ntracks()++;
          }
          vertices[thisVertex].isGood() = false;
        }
        bool closeCluster = (k + 1 == nCand);
        if (!closeCluster) {
          const int nextVertex = (int)vertices[k + 1].order();
          double combinedSigma = alpaka::math::sqrt(
              acc, sigma_z[thisVertex] * sigma_z[thisVertex] + sigma_z[nextVertex] * sigma_z[nextVertex]);
          double mergeThreshold = alpaka::math::max(acc, 2 * cParams.vertexSize, 3.0 * combinedSigma);
          closeCluster = alpaka::math::abs(acc, vertices[nextVertex].z() - vertices[thisVertex].z()) > mergeThreshold;
        }
        if (closeCluster) {
          vertices[head].isGood() = (vertices[head].ntracks() > 1);
          head = -1;
        }
      }
      // The merge (and the vertices without tracks) invalidated vertices, so the order list has to be compacted and the
      // vertex multiplicity set back. Done in one pass with a local counter, nV is written once at the end: shifting the
      // list in place while decrementing vertices[0].nV() inside the loop did not update nV (the loop bound kept its
      // initial value), which left nV about the number of removed vertices too large and the stale tail of the list
      // (copies of valid vertices) was then fitted again, giving duplicate vertices
      const int nOld = vertices[0].nV();
      int nKept = 0;
      for (int k = 0; k < nOld; k++) {
        const int thisVertex = vertices[k].order();
        if (thisVertex == emptyVertexOrder) {
          continue;  // purged, so it was bad
        }
        if (vertices[thisVertex].isGood()) {
          vertices[nKept].order() = thisVertex;
          nKept++;
        }
      }
      vertices[0].nV() = nKept;
    }
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ARBITRATOR
    if (once_per_block(acc)) {
      printf("[DAInBlocksClusterizerAlgoArbitrator::finalizeVertices()] Vertices are reordered, nV %i left\n",
             vertices[0].nV());
    }
#endif
    alpaka::syncBlockThreads(acc);
  }  //finalizeVertices

  class ArbitrateKernel {
  public:
    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC void operator()(const TAcc& acc,
                                  TrackForVertexDeviceCollection::View tracks,
                                  TrackForVertexScratchCollection::View scratch,
                                  VertexDeviceCollection::View vertices,
                                  DAInBlocksClusterParameters const cParams,
                                  int32_t nBlocks) const {
      // This has the core of the clusterization algorithm
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ARBITRATOR
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgoArbitrator::operator()] Start arbitration\n");
      }
#endif
      resortVerticesAndAssign(acc, tracks, scratch, vertices, cParams, nBlocks);
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ARBITRATOR
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgoArbitrator::operator()] Vertex reassignment finished\n");
        printf("[STAGE-CMP] after_resort nV=%d\n", vertices[0].nV());
      }
#endif
      alpaka::syncBlockThreads(acc);
      finalizeVertices(acc, tracks, scratch, vertices, cParams);  // In CUDA it used to be verticesAndClusterize
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ARBITRATOR
      if (once_per_block(acc)) {
        printf("[DAInBlocksClusterizerAlgoArbitrator::operator()] Vertices finalized\n");
        printf("[STAGE-CMP] arbitrate_final nV=%d\n", vertices[0].nV());
      }
#endif
      alpaka::syncBlockThreads(acc);
    }
  };  // class kernel

  void DAInBlocksClusterizerAlgo::arbitrate(Queue& queue,
                                            TrackForVertexDeviceCollection& deviceTrack,
                                            TrackForVertexScratchCollection& scratch,
                                            VertexDeviceCollection& deviceVertex,
                                            DAInBlocksClusterParameters const& cParams,
                                            int32_t nBlocks,
                                            int32_t blockSize) {
    const int blocks = 1;  // A single block: all the per-block vertex lists are merged into one collection here
    alpaka::exec<Acc1D>(queue,
                        make_workdiv<Acc1D>(blocks, blockSize),
                        ArbitrateKernel{},
                        deviceTrack.view(),
                        scratch.view(),
                        deviceVertex.view(),
                        cParams,
                        nBlocks);
  }  // arbitraterAlgo::arbitrate

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
