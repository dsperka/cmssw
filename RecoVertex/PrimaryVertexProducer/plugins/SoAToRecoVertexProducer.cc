#include "DataFormats/VertexSoA/interface/VertexHostCollection.h"
#include "DataFormats/VertexSoA/interface/TrackForVertexHostCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"

/**
   * This plugin takes the SoA portableVertex and converts them to reco::Vertex, for usage within other workflows
   * - consuming set of reco::Tracks and portablevertex SoA
   * - produces a host reco::vertexCollection
 */
class SoAToRecoVertexProducer : public edm::stream::EDProducer<> {
public:
  SoAToRecoVertexProducer(edm::ParameterSet const& config)
      : portableVertexToken_(consumes(config.getParameter<edm::InputTag>("soaVertex"))),
        // track_id() in the vertex SoA is an index into the per-block-duplicated tracksInBlocks
        // device collection (published by PrimaryVertexProducerPortable under the same module
        // label as soaVertex), not into srcTrack directly -- tracksInBlocks' own ttIndex()
        // column holds the true srcTrack index. See the comment in produce() below.
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
  void produce(edm::Event&, const edm::EventSetup&) override;
  const edm::EDGetTokenT<VertexHostCollection> portableVertexToken_;
  const edm::EDGetTokenT<TrackForVertexHostCollection> soaTracksInBlocksToken_;
  const edm::EDGetTokenT<reco::TrackCollection> recoTrackToken_;
  const edm::EDPutTokenT<reco::VertexCollection> recoVertexToken_;
};

void SoAToRecoVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Book inputs and space for outputs
  const VertexHostCollection& hostVertex = iEvent.get(portableVertexToken_);
  const VertexHostCollection::ConstView& hostVertexView = hostVertex.const_view();
  const TrackForVertexHostCollection& hostTracksInBlocks = iEvent.get(soaTracksInBlocksToken_);
  const TrackForVertexHostCollection::ConstView& hostTracksInBlocksView = hostTracksInBlocks.const_view();
  // Note that we need reco::Tracks for building the track Reference vector inside the reco::Vertex
  auto tracks = iEvent.getHandle(recoTrackToken_);

  auto result = std::make_unique<reco::VertexCollection>();

  // Do the conversion back to reco::Vertex
  reco::VertexCollection& vColl = (*result);
  // nV() good, compacted vertices are reached via the order() indirection (vertices[k].order()
  // for k=0..nV()-1 gives their RAW indices, scattered across the full buffer) -- NOT via raw
  // index k directly. Indexing hostVertexView[iV] with the loop counter itself (as this used to)
  // checks essentially arbitrary/unrelated raw slots for isGood(), silently dropping the large
  // majority of legitimately-finalized vertices from the output.
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
      // track_id() is an index into the tracksInBlocks device collection (deliberately -- see
      // the comment on PrimaryVertexProducerPortable::tracksInBlocksPutToken_), NOT an index
      // into srcTrack directly. Translate through ttIndex(), which tracksInBlocks carries
      // specifically to allow recovering the original srcTrack index. Building the TrackRef
      // directly from track_id() (as this used to) silently produced TrackRefs pointing at
      // unrelated/wrong (or out-of-range) entries of srcTrack.
      int blockLocalTrack = hostVertexView[iV].track_id()[iT];
      int new_itrack = hostTracksInBlocksView[blockLocalTrack].ttIndex();
      reco::TrackRef ref(tracks, new_itrack);
      newV.add(ref, hostVertexView[iV].track_weight()[iT]);
    }
  }
  // And finally put the collection in the event
  iEvent.put(std::move(result));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SoAToRecoVertexProducer);
