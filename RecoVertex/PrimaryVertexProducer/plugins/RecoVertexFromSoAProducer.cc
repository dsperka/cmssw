#include "DataFormats/VertexSoA/interface/VertexHostCollection.h"
#include "DataFormats/VertexSoA/interface/TrackForVertexHostCollection.h"
#include "DataFormats/VertexSoA/interface/VertexTrackHostCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexTimeAlgorithmFromTracksPID.h"
#include <memory>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"

#include <algorithm>

/**
   * This plugin takes the SoA portableVertex and converts them to reco::Vertex, for usage within other workflows
   * - consuming set of reco::Tracks and portablevertex SoA
   * - produces a host reco::vertexCollection
 */
class RecoVertexFromSoAProducer : public edm::stream::EDProducer<> {
public:
  RecoVertexFromSoAProducer(edm::ParameterSet const& config)
      : portableVertexToken_(consumes(config.getParameter<edm::InputTag>("soaVertex"))),
        soaTracksInBlocksToken_(consumes(config.getParameter<edm::InputTag>("soaVertex"))),
        soaVertexTracksToken_(consumes(config.getParameter<edm::InputTag>("soaVertex"))),
        recoTrackToken_(consumes(config.getParameter<edm::InputTag>("srcTrack"))),
        theTTBToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
        recoVertexToken_(produces<reco::VertexCollection>()) {
    // optional vertex time, computed as for the standard 3D producer (same algorithm, same inputs)
    const auto& timeConf = config.getParameter<edm::ParameterSet>("vertexTimeParameters");
    const std::string timeAlgo = timeConf.getParameter<std::string>("algorithm");
    if (timeAlgo == "fromTracksPID") {
      edm::ConsumesCollector collector = consumesCollector();
      timeEstimator_ = std::make_unique<VertexTimeAlgorithmFromTracksPID>(
          timeConf.getParameter<edm::ParameterSet>("fromTracksPID"), collector);
    } else if (!timeAlgo.empty()) {
      throw cms::Exception("Configuration")
          << "RecoVertexFromSoAProducer: unknown vertexTimeParameters.algorithm " << timeAlgo;
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("soaVertex");
    desc.add<edm::InputTag>("srcTrack");
    edm::ParameterSetDescription psdTime;
    psdTime.add<std::string>("algorithm", "");
    edm::ParameterSetDescription psdPID;
    VertexTimeAlgorithmFromTracksPID::fillPSetDescription(psdPID);
    psdTime.add<edm::ParameterSetDescription>("fromTracksPID", psdPID);
    desc.add<edm::ParameterSetDescription>("vertexTimeParameters", psdTime);

    descriptions.addWithDefaultLabel(desc);
  }

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  const edm::EDGetTokenT<::reco::VertexHostCollection> portableVertexToken_;
  const edm::EDGetTokenT<::reco::TrackForVertexHostCollection> soaTracksInBlocksToken_;
  const edm::EDGetTokenT<::reco::VertexTrackHostCollection> soaVertexTracksToken_;
  const edm::EDGetTokenT<reco::TrackCollection> recoTrackToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken_;
  const edm::EDPutTokenT<reco::VertexCollection> recoVertexToken_;
  std::unique_ptr<VertexTimeAlgorithmFromTracksPID> timeEstimator_;
};

void RecoVertexFromSoAProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Book inputs and space for outputs
  const ::reco::VertexHostCollection& hostVertex = iEvent.get(portableVertexToken_);
  const ::reco::VertexHostCollection::ConstView& hostVertexView = hostVertex.const_view();
  const ::reco::TrackForVertexHostCollection& hostTracksInBlocks = iEvent.get(soaTracksInBlocksToken_);
  const ::reco::TrackForVertexHostCollection::ConstView& hostTracksInBlocksView = hostTracksInBlocks.const_view();
  const ::reco::VertexTrackHostCollection::ConstView& hostVertexTracksView =
      iEvent.get(soaVertexTracksToken_).const_view();
  // Note that we need reco::Tracks for building the track Reference vector inside the reco::Vertex
  auto tracks = iEvent.getHandle(recoTrackToken_);

  const TransientTrackBuilder* ttBuilder = timeEstimator_ ? &iSetup.getData(theTTBToken_) : nullptr;
  if (timeEstimator_)
    timeEstimator_->setEvent(iEvent, iSetup);

  auto result = std::make_unique<reco::VertexCollection>();

  // Do the conversion back to reco::Vertex
  reco::VertexCollection& vColl = (*result);
  for (int iVO = 0; iVO < hostVertexView[0].nV(); iVO++) {
    int iV = hostVertexView[iVO].order();
    if (not(hostVertexView[iV].isGood()))
      continue;
    // Convert the SoA errors to a diagonal 3x3 matrix
    AlgebraicSymMatrix33 err;
    err[0][0] = hostVertexView[iV].errx();
    err[1][1] = hostVertexView[iV].erry();
    err[2][2] = hostVertexView[iV].errz();
    const int nTrk = hostVertexView[iV].ntracks();
    // reco::Track references used for building the vertex
    std::vector<reco::TrackRef> trkRefs;
    std::vector<float> trkWeights;
    trkRefs.reserve(nTrk);
    trkWeights.reserve(nTrk);
    for (int iT = 0; iT < nTrk; iT++) {
      const int iEntry = hostVertexView[iV].first_track() + iT;
      int blockLocalTrack = hostVertexTracksView[iEntry].track_id();
      int new_itrack = hostTracksInBlocksView[blockLocalTrack].ttIndex();
      trkRefs.emplace_back(tracks, new_itrack);
      trkWeights.push_back(hostVertexTracksView[iEntry].track_weight());
    }
    // Optional vertex time, from the same algorithm as the standard 3D producer
    double vtxTime = 0.;
    double vtxTimeError = 0.;
    bool hasTime = false;
    if (timeEstimator_) {
      std::vector<reco::TransientTrack> ttracks;
      ttracks.reserve(nTrk);
      for (const auto& ref : trkRefs)
        ttracks.push_back(ttBuilder->build(ref));
      TransientVertex tv(GlobalPoint(hostVertexView[iV].x(), hostVertexView[iV].y(), hostVertexView[iV].z()),
                         GlobalError(err(0, 0), 0, err(1, 1), 0, 0, err(2, 2)),
                         ttracks,
                         hostVertexView[iV].chi2(),
                         hostVertexView[iV].ndof());
      TransientVertex::TransientTrackToFloatMap wmap;
      for (int iT = 0; iT < nTrk; iT++)
        wmap[ttracks[iT]] = trkWeights[iT];
      tv.weightMap(wmap);
      float t = 0.f, terr = -1.f;
      if (timeEstimator_->vertexTime(t, terr, tv)) {
        vtxTime = t;
        vtxTimeError = terr;
        hasTime = true;
      }
    }
    // Then we can actually create the vertex
    if (hasTime) {
      reco::Vertex::Error4D err4;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j <= i; j++)
          err4(i, j) = err(i, j);
      err4(3, 3) = vtxTimeError * vtxTimeError;
      vColl.emplace_back(reco::Vertex::Point(hostVertexView[iV].x(), hostVertexView[iV].y(), hostVertexView[iV].z()),
                         err4,
                         vtxTime,
                         hostVertexView[iV].chi2(),
                         hostVertexView[iV].ndof(),
                         hostVertexView[iV].ntracks());
    } else {
      vColl.emplace_back(reco::Vertex::Point(hostVertexView[iV].x(), hostVertexView[iV].y(), hostVertexView[iV].z()),
                         err,
                         hostVertexView[iV].chi2(),
                         hostVertexView[iV].ndof(),
                         hostVertexView[iV].ntracks());
    }
    auto& newV = vColl.back();
    for (int iT = 0; iT < nTrk; iT++)
      newV.add(trkRefs[iT], trkWeights[iT]);
  }
  std::sort(vColl.begin(), vColl.end(), VertexHigherPtSquared());
  // And finally put the collection in the event
  iEvent.put(std::move(result));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoVertexFromSoAProducer);
