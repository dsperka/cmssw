#include "DataFormats/VertexSoA/interface/VertexHostCollection.h"
#include "DataFormats/VertexSoA/interface/TrackForVertexHostCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
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
class RecoVertexFromSoAProducer : public edm::global::EDProducer<> {
public:
  RecoVertexFromSoAProducer(edm::ParameterSet const& config)
      : portableVertexToken_(consumes(config.getParameter<edm::InputTag>("soaVertex"))),
        soaTracksInBlocksToken_(consumes(config.getParameter<edm::InputTag>("soaVertex"))),
        recoTrackToken_(consumes(config.getParameter<edm::InputTag>("srcTrack"))),
        recoVertexToken_(produces<reco::VertexCollection>()) {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("soaVertex");
    desc.add<edm::InputTag>("srcTrack");

    descriptions.addWithDefaultLabel(desc);
  }

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  const edm::EDGetTokenT<::reco::VertexHostCollection> portableVertexToken_;
  const edm::EDGetTokenT<::reco::TrackForVertexHostCollection> soaTracksInBlocksToken_;
  const edm::EDGetTokenT<reco::TrackCollection> recoTrackToken_;
  const edm::EDPutTokenT<reco::VertexCollection> recoVertexToken_;
};

void RecoVertexFromSoAProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  // Book inputs and space for outputs
  const ::reco::VertexHostCollection& hostVertex = iEvent.get(portableVertexToken_);
  const ::reco::VertexHostCollection::ConstView& hostVertexView = hostVertex.const_view();
  const ::reco::TrackForVertexHostCollection& hostTracksInBlocks = iEvent.get(soaTracksInBlocksToken_);
  const ::reco::TrackForVertexHostCollection::ConstView& hostTracksInBlocksView = hostTracksInBlocks.const_view();
  // Note that we need reco::Tracks for building the track Reference vector inside the reco::Vertex
  auto tracks = iEvent.getHandle(recoTrackToken_);

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
    // Then we can actually create the vertex
    vColl.emplace_back(reco::Vertex::Point(hostVertexView[iV].x(), hostVertexView[iV].y(), hostVertexView[iV].z()),
                       err,
                       hostVertexView[iV].chi2(),
                       hostVertexView[iV].ndof(),
                       hostVertexView[iV].ntracks());
    auto& newV = vColl.back();
    // Finally, add references to the reco::Track used for building it
    for (int iT = 0; iT < hostVertexView[iV].ntracks(); iT++) {
      int blockLocalTrack = hostVertexView[iV].track_id()[iT];
      int new_itrack = hostTracksInBlocksView[blockLocalTrack].ttIndex();
      reco::TrackRef ref(tracks, new_itrack);
      newV.add(ref, hostVertexView[iV].track_weight()[iT]);
    }
  }
  std::sort(vColl.begin(), vColl.end(), VertexHigherPtSquared());
  // And finally put the collection in the event
  iEvent.put(std::move(result));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoVertexFromSoAProducer);
