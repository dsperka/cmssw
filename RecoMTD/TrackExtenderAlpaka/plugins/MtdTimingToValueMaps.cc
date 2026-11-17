// Writes the 22 ValueMaps of the legacy TrackExtenderWithMTD (same instance names, keyed on the
// input tracks) from the MtdTiming table. Optionally (trackClones) also writes what the legacy module
// writes as its default products: one output track per input track (same parameters: there is no
// refit) carrying t0 and beta, with a TrackExtra whose hits are the track's tracker hits plus the
// chosen MTD hit(s). The MTD validation needs these hits.

#include <cmath>
#include <memory>
#include <vector>

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/MtdTimingSoA/interface/MtdTimingHost.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

class MtdTimingToValueMaps : public edm::global::EDProducer<> {
public:
  explicit MtdTimingToValueMaps(edm::ParameterSet const& cfg)
      : tracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("tracksSrc"))),
        timingToken_(consumes<mtdtiming::MtdTimingHost>(cfg.getParameter<edm::InputTag>("timingSrc"))),
        clones_(cfg.getParameter<bool>("trackClones")) {
    for (auto const& name : kFloatMaps)
      floatTokens_.push_back(produces<edm::ValueMap<float>>(name));
    npixBarrelToken_ = produces<edm::ValueMap<int>>("npixBarrel");
    mvaToken_ = produces<edm::ValueMap<float>>("mtdQualMVA");
    npixEndcapToken_ = produces<edm::ValueMap<int>>("npixEndcap");
    assocToken_ = produces<edm::ValueMap<int>>("generalTrackassoc");
    posToken_ = produces<edm::ValueMap<GlobalPoint>>("generalTrackmtdpos");
    if (clones_) {
      clonesToken_ = produces<reco::TrackCollection>();
      extrasToken_ = produces<reco::TrackExtraCollection>();
      hitsOutToken_ = produces<TrackingRecHitCollection>();
      hitsToken_ = consumes<MTDTrackingDetSetVector>(cfg.getParameter<edm::InputTag>("hitsSrc"));
      topoToken_ = esConsumes();
    }
  }

  void produce(edm::StreamID, edm::Event& ev, edm::EventSetup const& es) const override {
    auto const tracksH = ev.getHandle(tracksToken_);
    auto const& tracks = *tracksH;
    auto const& timing = ev.get(timingToken_);
    auto const v = timing.const_view();
    size_t const n = tracks.size();
    if (size_t(v.metadata().size()) != n)
      throw cms::Exception("MtdTimingToValueMaps") << "timing table has " << v.metadata().size() << " rows for " << n
                                                   << " tracks";

    std::vector<std::vector<float>> f(kFloatMaps.size(), std::vector<float>(n));
    std::vector<int> npixB(n), npixE(n), assoc(n);
    std::vector<float> mva(n);
    std::vector<GlobalPoint> pos(n);
    for (size_t i = 0; i < n; ++i) {
      auto const& trk = tracks[i];
      auto const r = v[i];
      bool const btl = r.btlHit() >= 0, etl = r.etlHit1() >= 0;
      float outermost;
      if (btl || etl)
        outermost = btl ? trk.outerRadius() : trk.outerZ();
      else
        outermost = std::abs(trk.eta()) < 1.5f ? trk.outerRadius() : trk.outerZ();
      float const vals[] = {r.btlChi2(),
                            r.etlChi2(),
                            r.btlTimeChi2(),
                            r.etlTimeChi2(),
                            outermost,
                            float(trk.p()),
                            r.beta(),
                            r.t0(),
                            r.sigmat0(),
                            r.pathLength(),
                            r.tmtd(),
                            r.sigmatmtd(),
                            r.tofPi(),
                            r.tofK(),
                            r.tofP(),
                            r.sigmaTofPi(),
                            r.sigmaTofK(),
                            r.sigmaTofP()};
      for (size_t k = 0; k < kFloatMaps.size(); ++k)
        f[k][i] = vals[k];
      npixB[i] = trk.hitPattern().numberOfValidPixelBarrelHits();
      npixE[i] = trk.hitPattern().numberOfValidPixelEndcapHits();
      assoc[i] = i;
      mva[i] = r.mva();
      pos[i] = GlobalPoint(r.mtdPosX(), r.mtdPosY(), r.mtdPosZ());
    }
    for (size_t k = 0; k < kFloatMaps.size(); ++k)
      fill(ev, tracksH, f[k], floatTokens_[k]);
    fill(ev, tracksH, npixB, npixBarrelToken_);
    fill(ev, tracksH, mva, mvaToken_);
    fill(ev, tracksH, npixE, npixEndcapToken_);
    fill(ev, tracksH, assoc, assocToken_);
    fill(ev, tracksH, pos, posToken_);

    if (clones_) {
      auto const& topo = es.getData(topoToken_);
      auto const& mtdHits = ev.get(hitsToken_).data();
      auto const hitsRefProd = ev.getRefBeforePut<TrackingRecHitCollection>(hitsOutToken_);
      auto const extrasRefProd = ev.getRefBeforePut<reco::TrackExtraCollection>(extrasToken_);
      auto out = std::make_unique<reco::TrackCollection>();
      auto extras = std::make_unique<reco::TrackExtraCollection>();
      auto outHits = std::make_unique<TrackingRecHitCollection>();
      out->reserve(n);
      for (size_t i = 0; i < n; ++i) {
        auto const& trk = tracks[i];
        auto const r = v[i];
        float const sigmat0 = r.sigmat0();
        out->emplace_back(trk.chi2(),
                          trk.ndof(),
                          trk.referencePoint(),
                          trk.momentum(),
                          trk.charge(),
                          trk.covariance(),
                          trk.algo(),
                          reco::TrackBase::undefQuality,
                          r.t0(),
                          r.beta(),
                          sigmat0 > 0.f ? sigmat0 * sigmat0 : -1.f,
                          r.timed() ? r.sigmaBeta() * r.sigmaBeta() : -1.f);
        auto& c = out->back();
        c.setQualityMask(trk.qualityMask());
        if (!trk.extra().isAvailable())
          continue;
        size_t const first = outHits->size();
        for (auto const& hit : trk.recHits())
          outHits->push_back(hit->clone());
        for (int32_t idx : {r.btlHit(), r.etlHit1(), r.etlHit2()})
          if (idx >= 0 && size_t(idx) < mtdHits.size())
            outHits->push_back(mtdHits[idx].clone());
        auto const& ox = *trk.extra();
        extras->emplace_back(ox.outerPosition(),
                             ox.outerMomentum(),
                             ox.outerOk(),
                             ox.innerPosition(),
                             ox.innerMomentum(),
                             ox.innerOk(),
                             ox.outerStateCovariance(),
                             ox.outerDetId(),
                             ox.innerStateCovariance(),
                             ox.innerDetId(),
                             ox.seedDirection(),
                             ox.seedRef());
        extras->back().setHits(hitsRefProd, first, outHits->size() - first);
        c.setExtra(reco::TrackExtraRef(extrasRefProd, extras->size() - 1));
        for (size_t h = first; h < outHits->size(); ++h)
          c.appendHitPattern((*outHits)[h], topo);
      }
      ev.put(clonesToken_, std::move(out));
      ev.put(extrasToken_, std::move(extras));
      ev.put(hitsOutToken_, std::move(outHits));
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracksSrc", edm::InputTag("generalTracks"));
    desc.add<edm::InputTag>("timingSrc", edm::InputTag("mtdTimingProducer"));
    desc.add<bool>("trackClones", true);
    desc.add<edm::InputTag>("hitsSrc", edm::InputTag("mtdTrackingRecHits"));
    descriptions.addWithDefaultLabel(desc);
  }

private:
  template <class T>
  static void fill(edm::Event& ev,
                   edm::Handle<reco::TrackCollection> const& h,
                   std::vector<T> const& vec,
                   edm::EDPutToken const& token) {
    auto out = std::make_unique<edm::ValueMap<T>>();
    typename edm::ValueMap<T>::Filler filler(*out);
    filler.insert(h, vec.begin(), vec.end());
    filler.fill();
    ev.put(token, std::move(out));
  }

  inline static const std::vector<std::string> kFloatMaps = {"btlMatchChi2",
                                                             "etlMatchChi2",
                                                             "btlMatchTimeChi2",
                                                             "etlMatchTimeChi2",
                                                             "generalTrackOutermostHitPosition",
                                                             "generalTrackp",
                                                             "generalTrackBeta",
                                                             "generalTrackt0",
                                                             "generalTracksigmat0",
                                                             "generalTrackPathLength",
                                                             "generalTracktmtd",
                                                             "generalTracksigmatmtd",
                                                             "generalTrackTofPi",
                                                             "generalTrackTofK",
                                                             "generalTrackTofP",
                                                             "generalTrackSigmaTofPi",
                                                             "generalTrackSigmaTofK",
                                                             "generalTrackSigmaTofP"};

  edm::EDGetTokenT<reco::TrackCollection> const tracksToken_;
  edm::EDGetTokenT<mtdtiming::MtdTimingHost> const timingToken_;
  bool const clones_;
  std::vector<edm::EDPutToken> floatTokens_;
  edm::EDPutToken mvaToken_;
  edm::EDPutToken npixBarrelToken_, npixEndcapToken_, assocToken_, posToken_, clonesToken_, extrasToken_, hitsOutToken_;
  edm::EDGetTokenT<MTDTrackingDetSetVector> hitsToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
};

DEFINE_FWK_MODULE(MtdTimingToValueMaps);
