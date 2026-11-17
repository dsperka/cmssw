#ifndef RecoMTD_TrackExtenderAlpaka_plugins_alpaka_TofPidAlgo_h
#define RecoMTD_TrackExtenderAlpaka_plugins_alpaka_TofPidAlgo_h

#include "DataFormats/MtdTimingSoA/interface/MtdEventInputSoA.h"
#include "DataFormats/MtdTimingSoA/interface/MtdTimingSoA.h"
#include "DataFormats/MtdTimingSoA/interface/alpaka/TofPidDevice.h"
#include "DataFormats/VertexSoA/interface/TrackForVertexSoA.h"
#include "DataFormats/VertexSoA/interface/VertexSoA.h"
#include "DataFormats/VertexSoA/interface/VertexTrackSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming {

  // same names and meaning as the TOFPIDProducer parameters
  struct TofPidParameters {
    double vtxMaxSigmaT;
    double maxDz;
    double maxDtSignificance;
    double minProbHeavy;
    double fixedT0Error;
    double probPion, probKaon, probProton;
    bool vertexReassignment;
    bool MVASel;
    double minTrackTimeQuality;
  };

  void runTofPid(Queue& queue,
                 TofPidParameters const& params,
                 ::mtdtiming::MtdTrackInSoA::ConstView tracks,
                 ::mtdtiming::MtdTimingSoA::ConstView timing,
                 ::reco::VertexSoA::ConstView vertices,
                 ::reco::VertexTrackSoA::ConstView vertexTracks,
                 ::reco::TrackForVertexSoA::ConstView tracksInBlocks,
                 TofPidDevice& out);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming

#endif
