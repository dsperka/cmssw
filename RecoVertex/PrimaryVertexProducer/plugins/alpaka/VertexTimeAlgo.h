#ifndef RecoVertex_PrimaryVertexProducer_plugins_alpaka_VertexTimeAlgo_h
#define RecoVertex_PrimaryVertexProducer_plugins_alpaka_VertexTimeAlgo_h

// Vertex time from the MTD times of its tracks, on the device: the same algorithm as
// VertexTimeAlgorithmFromTracksPID (pion/kaon/proton hypotheses, annealing in 1/T).

#include "DataFormats/MtdTimingSoA/interface/MtdTimingSoA.h"
#include "DataFormats/VertexSoA/interface/alpaka/TrackForVertexDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexTrackDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  struct VertexTimeParameters {
    double minTrackVtxWeight;    // 0.5
    double minTrackTimeQuality;  // 0.8
    double probPion, probKaon, probProton;
    double Tstart;         // 256
    double coolingFactor;  // 0.5
    bool useMVAVtxTime;    // true
  };

  // Fills t and errt of every good vertex; errt = -1 if the vertex gets no time.
  void vertexTimes(Queue& queue,
                   VertexTimeParameters const& params,
                   ::reco::TrackForVertexSoA::ConstView tracksInBlocks,
                   ::reco::VertexTrackSoA::ConstView vertexTracks,
                   ::mtdtiming::MtdTimingSoA::ConstView timing,
                   ::reco::VertexSoA::View vertices);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
