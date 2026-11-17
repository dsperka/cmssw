#include <algorithm>
#include <type_traits>

#include "DataFormats/VertexSoA/interface/alpaka/VertexDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/TrackForVertexDeviceCollection.h"
#include "DataFormats/Portable/interface/PortableCollection.h"
#include "DataFormats/VertexSoA/interface/VertexHostCollection.h"
#include "DataFormats/VertexSoA/interface/TrackForVertexHostCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackForVertexFiller.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  /**
   * This:
   * - consumes set of reco::Tracks and reco::BeamSpot
   * - converts the reco::Tracks to a Alpaka-friendly dataformat ::reco::TrackForVertexHostCollection
   * - puts the Alpaka dataformat in the device for later consumption
   */
  class PortableTrackSoAProducer : public global::EDProducer<> {
  public:
    PortableTrackSoAProducer(edm::ParameterSet const& config)
        : global::EDProducer<>(config),
          theTTBToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
          trackToken_(consumes<reco::TrackCollection>(config.getParameter<edm::InputTag>("TrackLabel"))),
          beamSpotToken_(consumes<reco::BeamSpot>(config.getParameter<edm::InputTag>("BeamSpotLabel"))),
          fParams(config.getParameter<edm::ParameterSet>("TkFilterParameters")) {
      hostPutToken_ = produces();
    }

    void produce(edm::StreamID sid, device::Event& iEvent, device::EventSetup const& iSetup) const override {
      // Get input collections from event
      auto tracks = iEvent.getHandle(trackToken_);
      const reco::BeamSpot& beamSpot = iEvent.get(beamSpotToken_);
      int32_t tsize_ = tracks->size();

      ::reco::TrackForVertexHostCollection hostTracks(iEvent.queue(), tsize_);
      pvalpaka::fillTracksForVertexing(hostTracks.view(), tracks, beamSpot, iSetup.getData(theTTBToken_), fParams);
      // Put the host collection into the event: consumers that ask for the device collection get it copied
      // implicitly by the framework (and no copy at all on the CPU backend)
      iEvent.emplace(hostPutToken_, std::move(hostTracks));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("TrackLabel");
      desc.add<edm::InputTag>("BeamSpotLabel");
      auto psd0 = pvalpaka::TrackFilterParametersForVertexing::description();
      desc.add<edm::ParameterSetDescription>("TkFilterParameters", psd0);
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken_;
    const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
    const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    edm::EDPutTokenT<::reco::TrackForVertexHostCollection> hostPutToken_;
    const pvalpaka::TrackFilterParametersForVertexing fParams;
  };  //PortableTrackSoAProducer declaration

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(PortableTrackSoAProducer);
