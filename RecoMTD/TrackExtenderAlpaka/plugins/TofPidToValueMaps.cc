// Host side of the device TOF particle identification: the TofPid table to the ValueMaps of
// TOFPIDProducer (same product names), so that downstream modules read it unchanged.

#include <array>
#include <vector>

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/MtdTimingSoA/interface/TofPidHost.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

class TofPidToValueMaps : public edm::global::EDProducer<> {
public:
  explicit TofPidToValueMaps(edm::ParameterSet const& cfg)
      : tracksToken_{consumes(cfg.getParameter<edm::InputTag>("tracksSrc"))},
        pidToken_{consumes(cfg.getParameter<edm::InputTag>("tofPidSrc"))} {
    for (auto const* name : kNames)
      tokens_.push_back(produces<edm::ValueMap<float>>(name));
  }

  void produce(edm::StreamID, edm::Event& ev, edm::EventSetup const&) const override {
    auto const tracksH = ev.getHandle(tracksToken_);
    auto const v = ev.get(pidToken_).const_view();
    size_t const n = tracksH->size();
    if (size_t(v.metadata().size()) != n)
      throw cms::Exception("TofPidToValueMaps") << "TofPid table has " << v.metadata().size() << " rows for " << n
                                                << " tracks";
    std::vector<std::vector<float>> f(kNames.size(), std::vector<float>(n));
    for (size_t i = 0; i < n; ++i) {
      auto const r = v[i];
      float const vals[] = {r.t0(), r.sigmat0(), r.t0safe(), r.sigmat0safe(), r.probPi(), r.probK(), r.probP()};
      for (size_t k = 0; k < kNames.size(); ++k)
        f[k][i] = vals[k];
    }
    for (size_t k = 0; k < kNames.size(); ++k) {
      auto out = std::make_unique<edm::ValueMap<float>>();
      edm::ValueMap<float>::Filler filler(*out);
      filler.insert(tracksH, f[k].begin(), f[k].end());
      filler.fill();
      ev.put(tokens_[k], std::move(out));
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracksSrc", edm::InputTag("generalTracks"));
    desc.add<edm::InputTag>("tofPidSrc", edm::InputTag("tofPidProducer"));
    descriptions.addWithDefaultLabel(desc);
  }

private:
  static constexpr std::array<char const*, 7> kNames = {"t0", "sigmat0", "t0safe", "sigmat0safe", "probPi", "probK", "probP"};
  edm::EDGetTokenT<reco::TrackCollection> const tracksToken_;
  edm::EDGetTokenT<mtdtiming::TofPidHost> const pidToken_;
  std::vector<edm::EDPutTokenT<edm::ValueMap<float>>> tokens_;
};

DEFINE_FWK_MODULE(TofPidToValueMaps);
