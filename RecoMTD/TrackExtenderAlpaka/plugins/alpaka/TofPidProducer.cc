// TOF particle identification (TOFPIDProducer) on the device, from the MTD timing of the tracks and the
// device vertices (PrimaryVertexProducerPortable with timingSrc set, so that they have a time).
// Output: one TofPid row per input track. No host-device copy is made here.

#include "DataFormats/MtdTimingSoA/interface/alpaka/MtdEventInputDevice.h"
#include "DataFormats/MtdTimingSoA/interface/alpaka/MtdTimingDevice.h"
#include "DataFormats/MtdTimingSoA/interface/alpaka/TofPidDevice.h"
#include "DataFormats/VertexSoA/interface/alpaka/TrackForVertexDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexTrackDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include "TofPidAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class TofPidProducer : public global::EDProducer<> {
  public:
    explicit TofPidProducer(edm::ParameterSet const& cfg)
        : EDProducer<>(cfg),
          inputToken_{consumes(cfg.getParameter<edm::InputTag>("inputSrc"))},
          timingToken_{consumes(cfg.getParameter<edm::InputTag>("timingSrc"))},
          vertexToken_{consumes(cfg.getParameter<edm::InputTag>("vtxsSrc"))},
          vertexTracksToken_{consumes(cfg.getParameter<edm::InputTag>("vtxsSrc"))},
          // the tracks in blocks are shared by the "" and "WithBS" fits: no instance name
          tracksInBlocksToken_{consumes(edm::InputTag(cfg.getParameter<edm::InputTag>("vtxsSrc").label()))},
          outToken_{produces()},
          params_{.vtxMaxSigmaT = cfg.getParameter<double>("vtxMaxSigmaT"),
                  .maxDz = cfg.getParameter<double>("maxDz"),
                  .maxDtSignificance = cfg.getParameter<double>("maxDtSignificance"),
                  .minProbHeavy = cfg.getParameter<double>("minProbHeavy"),
                  .fixedT0Error = cfg.getParameter<double>("fixedT0Error"),
                  .probPion = cfg.getParameter<double>("probPion"),
                  .probKaon = cfg.getParameter<double>("probKaon"),
                  .probProton = cfg.getParameter<double>("probProton"),
                  .vertexReassignment = cfg.getParameter<bool>("vertexReassignment"),
                  .MVASel = cfg.getParameter<bool>("MVASel"),
                  .minTrackTimeQuality = cfg.getParameter<double>("minTrackTimeQuality")} {}

    void produce(edm::StreamID, device::Event& ev, device::EventSetup const&) const override {
      auto const& input = ev.get(inputToken_);
      auto const& timing = ev.get(timingToken_);
      auto const tracks = input.const_view().tracks();
      mtdtiming::TofPidDevice out(ev.queue(), tracks.metadata().size());
      mtdtiming::runTofPid(ev.queue(),
                           params_,
                           tracks,
                           timing.const_view(),
                           ev.get(vertexToken_).const_view(),
                           ev.get(vertexTracksToken_).const_view(),
                           ev.get(tracksInBlocksToken_).const_view(),
                           out);
      ev.emplace(outToken_, std::move(out));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("inputSrc", edm::InputTag("mtdEventInput"));
      desc.add<edm::InputTag>("timingSrc", edm::InputTag("mtdTimingProducer"));
      desc.add<edm::InputTag>("vtxsSrc", edm::InputTag("pvVertexSoA"));
      // same names and defaults as TOFPIDProducer
      desc.add<double>("vtxMaxSigmaT", 0.025);
      desc.add<double>("maxDz", 0.1);
      desc.add<double>("maxDtSignificance", 5.0);
      desc.add<double>("minProbHeavy", 0.75);
      desc.add<double>("fixedT0Error", 0.);
      desc.add<double>("probPion", 1.);
      desc.add<double>("probKaon", 1.);
      desc.add<double>("probProton", 1.);
      desc.add<bool>("vertexReassignment", true);
      desc.add<bool>("MVASel", false);
      desc.add<double>("minTrackTimeQuality", 0.8);
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    device::EDGetToken<mtdtiming::MtdEventInputDevice> const inputToken_;
    device::EDGetToken<mtdtiming::MtdTimingDevice> const timingToken_;
    device::EDGetToken<VertexDeviceCollection> const vertexToken_;
    device::EDGetToken<VertexTrackDeviceCollection> const vertexTracksToken_;
    device::EDGetToken<TrackForVertexDeviceCollection> const tracksInBlocksToken_;
    device::EDPutToken<mtdtiming::TofPidDevice> const outToken_;
    mtdtiming::TofPidParameters const params_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(TofPidProducer);
