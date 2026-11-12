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

//#define DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ALPAKA_PORTABLETRACKSOAPRODUCER 0

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  /**
   * This:
   * - consumes set of reco::Tracks and reco::BeamSpot
   * - converts the reco::Tracks to a Alpaka-friendly dataformat ::reco::TrackForVertexHostCollection
   * - puts the Alpaka dataformat in the device for later consumption
   */
  struct TrackFilterParametersForVertexing {
    explicit TrackFilterParametersForVertexing(edm::ParameterSet const& pset)
        : maxSignificance(pset.getParameter<double>("maxD0Significance")),
          maxdxyError(pset.getParameter<double>("maxD0Error")),
          maxdzError(pset.getParameter<double>("maxDzError")),
          minpAtIP(pset.getParameter<double>("minPt")),
          maxetaAtIP(pset.getParameter<double>("maxEta")),
          maxchi2(pset.getParameter<double>("maxNormalizedChi2")),
          minpixelHits(pset.getParameter<int>("minPixelLayersWithHits")),
          mintrackerHits(pset.getParameter<int>("minSiliconLayersWithHits")),
          trackQuality(reco::TrackBase::undefQuality),
          vertexSize(pset.getParameter<double>("vertexSize")),
          d0CutOff(pset.getParameter<double>("d0CutOff")) {
      std::string const qualityClass = pset.getParameter<std::string>("trackQuality");
      if (qualityClass != "any" && qualityClass != "Any" && qualityClass != "ANY" && !(qualityClass.empty()))
        trackQuality = reco::TrackBase::qualityByName(qualityClass);
    }

    // Configurable filter parameters for the tracks
    double maxSignificance;
    double maxdxyError;
    double maxdzError;
    double minpAtIP;
    double maxetaAtIP;
    double maxchi2;
    int minpixelHits;
    int mintrackerHits;
    reco::TrackBase::TrackQuality trackQuality;
    double vertexSize;
    double d0CutOff;
  };

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

      // Host collections
      ::reco::TrackForVertexHostCollection hostTracks(iEvent.queue(), tsize_);
      auto& tview = hostTracks.view();

      // Fill Host collections with input, first initialize globals
      tview.totweight() = 0;
      tview.nT() = 0;

      // Build transient tracks
      const auto& theB = iSetup.getData(theTTBToken_);
      std::vector<reco::TransientTrack> t_tks = theB.build(tracks, beamSpot);

      // Sort by z at the beam line. stateAtBeamLine() returns by value, so evaluate z once per track rather than
      // twice per comparison. The index into the reco::Track collection travels along, to later redo the
      // conversion back to reco::Vertex.
      using ZType = std::decay_t<decltype(t_tks.front().stateAtBeamLine().trackStateAtPCA().position().z())>;
      std::vector<std::pair<ZType, int32_t>> sortedZAndIndex;
      sortedZAndIndex.reserve(tsize_);
      for (int32_t idx = 0; idx < tsize_; idx++) {
        sortedZAndIndex.emplace_back(t_tks[idx].stateAtBeamLine().trackStateAtPCA().position().z(), idx);
      }
      std::sort(sortedZAndIndex.begin(), sortedZAndIndex.end(), [](auto const& a, auto const& b) {
        return a.first < b.first;
      });
      // This will keep track of how many we actually copy to device, only those that pass filter
      int32_t nTrueTracks = 0;
      for (int32_t idx = 0; idx < tsize_; idx++) {
        // Fill up the the Track SoA, weight doubles up as an isGood flag, as we compute it only for good tracks
        const int32_t trackIdx = sortedZAndIndex[idx].second;
        double weight = convertTrack(tview[nTrueTracks], t_tks[trackIdx], beamSpot, trackIdx, nTrueTracks);
        if (weight > 0) {
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ALPAKA_PORTABLETRACKSOAPRODUCER
          printf("[PortableTrackSoAProducer::produce()] Add track at z=%1.5f \n", sortedZAndIndex[idx].first);
#endif
          nTrueTracks += 1;
          tview.nT() += 1;
          tview.totweight() += weight;
        }
      }
#ifdef DEBUG_RECOVERTEX_PRIMARYVERTEXPRODUCER_ALPAKA_PORTABLETRACKSOAPRODUCER
      printf("[PortableTrackSoAProducer::produce()] From %i tracks, %i pass filters\n",
             (int32_t)tracks->size(),
             nTrueTracks);
#endif
      // Put the host collection into the event: consumers that ask for the device collection get it copied
      // implicitly by the framework (and no copy at all on the CPU backend)
      iEvent.emplace(hostPutToken_, std::move(hostTracks));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("TrackLabel");
      desc.add<edm::InputTag>("BeamSpotLabel");
      edm::ParameterSetDescription psd0;
      psd0.add<double>("maxNormalizedChi2", 10.0);
      psd0.add<double>("minPt", 0.0);
      psd0.add<std::string>("algorithm", "filter");
      psd0.add<double>("maxEta", 2.4);
      psd0.add<double>("maxD0Significance", 4.0);
      psd0.add<double>("maxD0Error", 1.0);
      psd0.add<double>("maxDzError", 1.0);
      psd0.add<std::string>("trackQuality", "any");
      psd0.add<int>("minPixelLayersWithHits", 2);
      psd0.add<int>("minSiliconLayersWithHits", 5);
      psd0.add<double>("vertexSize", 0.006);
      psd0.add<double>("d0CutOff", 0.10);
      desc.add<edm::ParameterSetDescription>("TkFilterParameters", psd0);
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken_;
    const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
    const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    edm::EDPutTokenT<::reco::TrackForVertexHostCollection> hostPutToken_;
    double convertTrack(::reco::TrackForVertexHostCollection::View::element out,
                        const reco::TransientTrack& in,
                        const reco::BeamSpot& bs,
                        int32_t idx,
                        int32_t order) const;
    const TrackFilterParametersForVertexing fParams;
  };  //PortableTrackSoAProducer declaration

  double PortableTrackSoAProducer::convertTrack(::reco::TrackForVertexHostCollection::View::element out,
                                                const reco::TransientTrack& in,
                                                const reco::BeamSpot& bs,
                                                int32_t idx,
                                                int32_t order) const {
    double weight = -1;
    // stateAtBeamLine() returns by value, and so do the quantities derived from it: evaluate each of them once
    const auto stateAtBeamLine = in.stateAtBeamLine();
    const auto transverseIP = stateAtBeamLine.transverseImpactParameter();
    // First check if it passes filters
    if ((transverseIP.significance() < fParams.maxSignificance) && (transverseIP.error() < fParams.maxdxyError) &&
        (in.track().dzError() < fParams.maxdzError) &&
        (in.impactPointState().globalMomentum().transverse() > fParams.minpAtIP) &&
        (std::fabs(in.impactPointState().globalMomentum().eta()) < fParams.maxetaAtIP) &&
        (in.normalizedChi2() < fParams.maxchi2) &&
        (in.hitPattern().pixelLayersWithMeasurement() >= fParams.minpixelHits) &&
        (in.hitPattern().trackerLayersWithMeasurement() >= fParams.mintrackerHits) &&
        (in.track().quality(fParams.trackQuality) || (fParams.trackQuality == reco::TrackBase::undefQuality))) {
      weight = 1.;
      if (fParams.d0CutOff > 0) {
        // significance is measured in the transverse plane
        double significance = transverseIP.value() / transverseIP.error();
        // weight is based on transverse displacement of the track
        weight = 1. / (1. + exp(std::pow(significance, 2) - std::pow(fParams.d0CutOff, 2)));
      }
      const auto trackStateAtPCA = stateAtBeamLine.trackStateAtPCA();
      const auto position = trackStateAtPCA.position();
      const auto momentum = trackStateAtPCA.momentum();
      // Just fill up variables
      out.x() = position.x();
      out.y() = position.y();
      out.z() = position.z();
      out.px() = momentum.x();
      out.py() = momentum.y();
      out.pz() = momentum.z();
      out.weight() = weight;
      // The original index in the reco::Track collection so we can go back to it eventually
      out.ttIndex() = idx;
      out.dz2() = std::pow(in.track().dzError(), 2);
      // Modified dz2 to account correlations and vertex size for clusterizer
      // dz^2 + (bs*pt)^2*pz^2/pt^2 + vertexSize^2
      double oneoverdz2 = (out.dz2()) +
                          ((std::pow(bs.BeamWidthX() * out.px(), 2)) + (std::pow(bs.BeamWidthY() * out.py(), 2))) *
                              std::pow(out.pz(), 2) / std::pow(momentum.perp2(), 2) +
                          std::pow(fParams.vertexSize, 2);
      oneoverdz2 = 1. / oneoverdz2;
      out.oneoverdz2() = oneoverdz2;
      out.dxy2AtIP() = std::pow(in.track().dxyError(), 2);
      out.dxy2() = std::pow(transverseIP.error(), 2);
      out.order() = order;
      out.isGood() = true;  // if we are here, we are to keep this track
    }
    return weight;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(PortableTrackSoAProducer);
