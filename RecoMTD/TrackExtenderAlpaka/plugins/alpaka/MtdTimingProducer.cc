// MTD track timing, Alpaka version (no refit).
//
// Reads the packed event input (MtdEventInputProducer) and the BTL geometry (MtdBtlGeometryESProducer),
// both already on the device; the device part (MtdTimingAlgo.dev.cc) finds the best hit per layer and
// computes the timing. Output: one MtdTiming row per input track. No host-device copy is made here:
// the few numbers the kernels take as arguments are read from the host copies of the same products.

#include "CommonTools/MVAUtils/interface/GBRForestTools.h"
#include "CondFormats/GBRForest/interface/GBRForest.h"
#include "DataFormats/MtdTimingSoA/interface/MtdEventInputHost.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "HeterogeneousCore/AlpakaCore/interface/MoveToDeviceCache.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "DataFormats/MtdTimingSoA/interface/alpaka/MtdEventInputDevice.h"
#include "DataFormats/MtdTimingSoA/interface/alpaka/MtdTimingDevice.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"
#include "RecoMTD/TrackExtenderAlpaka/interface/MtdBtlGeometryHost.h"
#include "RecoMTD/TrackExtenderAlpaka/interface/alpaka/MtdBtlGeometryDevice.h"

#include "MtdTimingAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace ::mtdtiming;

  class MtdTimingProducer : public global::EDProducer<> {
  public:
    explicit MtdTimingProducer(edm::ParameterSet const& cfg)
        : EDProducer<>(cfg),
          inputToken_{consumes(cfg.getParameter<edm::InputTag>("inputSrc"))},
          inputHostToken_{consumes(cfg.getParameter<edm::InputTag>("inputSrc"))},
          geomToken_{esConsumes()},
          geomHostToken_{esConsumes()},
          outToken_{produces()},
          bdt_{loadBdt(cfg.getParameter<edm::FileInPath>("qualityBDT_weights_file"))} {
      auto const& m = cfg.getParameterSet("MTDHitMatcher");
      params_.btlChi2Cut = m.getParameter<double>("btlChi2Cut");
      params_.btlTimeChi2Cut = m.getParameter<double>("btlTimeChi2Cut");
      params_.etlChi2Cut = m.getParameter<double>("etlChi2Cut");
      params_.etlTimeChi2Cut = m.getParameter<double>("etlTimeChi2Cut");
      params_.estimatorMaxChi2 = m.getParameter<double>("estimatorMaxChi2");
      params_.nSigmaWindow = m.getParameter<double>("estimatorMaxNSigma");
      params_.bsTimeSpread = m.getParameter<double>("bsTimeSpread");
    }

    void produce(edm::StreamID, device::Event& ev, device::EventSetup const& es) const override {
      auto const& input = ev.get(inputToken_);
      auto const hostInput = ev.get(inputHostToken_).const_view();
      auto const& geometry = es.getData(geomToken_);
      auto const btl = es.getData(geomHostToken_).const_view().lookup();
      int32_t const nTracks = hostInput.tracks().metadata().size();

      // per-event numbers passed to the kernels by value
      MatchParams params = params_;
      params.btlHalfLength = btl.halfLength();
      params.btlHalfThickness = btl.halfThickness();
      params.btlNRods = es.getData(geomHostToken_).const_view().rods().metadata().size();
      params.btlRodPhiStep = btl.rodPhiStep();
      params.btlRodPhiOffset = btl.rodPhiOffset();
      params.btlRodsOverlap = btl.rodsOverlap();
      params.btlRodsPeriodic = btl.rodsPeriodic();
      auto const layers = hostInput.layers();
      for (int l = 0; l < kNLayers; ++l) {
        params.layerRef[l] = layers[l].ref();
        params.layerSpread[l] = layers[l].spread();
        params.maxHitErr[l] = layers[l].maxHitErr();
      }

      auto& queue = ev.queue();
      mtdtiming::RecordDevice records(queue, nTracks * kNLayers);
      mtdtiming::MtdTimingDevice out(queue, nTracks);
      mtdtiming::runMtdTiming(
          queue, params, input.const_view(), geometry.const_view(), bdt_.get(queue).const_view(), records, out);
      ev.emplace(outToken_, std::move(out));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("inputSrc", edm::InputTag("mtdEventInput"));
      edm::ParameterSetDescription m;  // same names and defaults as the legacy MTDHitMatcher PSet
      m.add<double>("estimatorMaxChi2", 500.);
      m.add<double>("estimatorMaxNSigma", 10.);
      m.add<double>("btlChi2Cut", 50.);
      m.add<double>("btlTimeChi2Cut", 10.);
      m.add<double>("etlChi2Cut", 50.);
      m.add<double>("etlTimeChi2Cut", 10.);
      m.add<double>("bsTimeSpread", 0.2);
      desc.add<edm::ParameterSetDescription>("MTDHitMatcher", m);
      desc.add<edm::FileInPath>("qualityBDT_weights_file",
                                edm::FileInPath("RecoMTD/TimingIDTools/data/BDT_nvars_17_d7.xml"))
          ->setComment("MTD track-quality BDT, as mtdTrackQualityMVA");
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    using BdtHost = PortableHostCollection<MtdBdtSoA>;

    // the GBRForest of the weight file as flat tables (same trees, same order)
    static BdtHost loadBdt(edm::FileInPath const& file) {
      auto const forest = createGBRForest(file);
      int32_t nNodes = 0, nResponses = 0;
      for (auto const& tree : forest->Trees()) {
        nNodes += tree.CutIndices().size();
        nResponses += tree.Responses().size();
      }
      int32_t const nTrees = forest->Trees().size();
      BdtHost bdt(cms::alpakatools::host(), nTrees, nNodes, nResponses);
      auto trees = bdt.view().trees();
      auto nodes = bdt.view().nodes();
      auto responses = bdt.view().responses();
      // GBRForest has no getter for its initial response (a protected member)
      struct Access : GBRForest {
        static double initialResponse(GBRForest const& f) { return f.*(&Access::fInitialResponse); }
      };
      trees.initialResponse() = Access::initialResponse(*forest);
      int32_t n = 0, r = 0;
      for (int32_t t = 0; t < nTrees; ++t) {
        auto const& tree = forest->Trees()[t];
        trees[t].firstNode() = n;
        trees[t].firstResponse() = r;
        for (size_t k = 0; k < tree.CutIndices().size(); ++k, ++n) {
          nodes[n].cutIndex() = tree.CutIndices()[k];
          nodes[n].cutVal() = tree.CutVals()[k];
          nodes[n].left() = tree.LeftIndices()[k];
          nodes[n].right() = tree.RightIndices()[k];
        }
        for (float const resp : tree.Responses())
          responses[r++].response() = resp;
      }
      return bdt;
    }

    device::EDGetToken<mtdtiming::MtdEventInputDevice> const inputToken_;
    edm::EDGetTokenT<MtdEventInputHost> const inputHostToken_;
    device::ESGetToken<mtdtiming::MtdBtlGeometryDevice, MTDRecoGeometryRecord> const geomToken_;
    edm::ESGetToken<MtdBtlGeometryHost, MTDRecoGeometryRecord> const geomHostToken_;
    device::EDPutToken<mtdtiming::MtdTimingDevice> const outToken_;
    MatchParams params_;
    cms::alpakatools::MoveToDeviceCache<Device, BdtHost> const bdt_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(MtdTimingProducer);
