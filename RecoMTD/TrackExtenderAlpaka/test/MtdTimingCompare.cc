// Validation dump: for every track, the legacy and the new timing maps side by side plus MC truth
// (best matching TrackingParticle: pdgId, production-vertex time). One CSV line per track.

#include <cstdio>
#include <string>
#include <vector>

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/Associations/interface/TrackAssociation.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"

namespace {
  const std::vector<std::string> kF = {"btlMatchChi2",
                                       "etlMatchChi2",
                                       "btlMatchTimeChi2",
                                       "etlMatchTimeChi2",
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
                                       "generalTrackSigmaTofP",
                                       "generalTrackOutermostHitPosition"};
  const std::vector<std::string> kShort = {"btlChi2", "etlChi2", "btlTChi2", "etlTChi2", "p",     "beta",
                                           "t0",      "sigmat0", "path",     "tmtd",     "stmtd", "tofpi",
                                           "tofk",    "tofp",    "stofpi",   "stofp",    "outer"};
}  // namespace

class MtdTimingCompare : public edm::one::EDAnalyzer<> {
public:
  explicit MtdTimingCompare(edm::ParameterSet const& cfg)
      : tracksToken_(consumes(cfg.getParameter<edm::InputTag>("tracks"))),
        r2sToken_(consumes(cfg.getParameter<edm::InputTag>("recoToSim"))),
        out_(std::fopen(cfg.getParameter<std::string>("fileName").c_str(), "w")) {
    for (auto const& label : {cfg.getParameter<std::string>("legacy"), cfg.getParameter<std::string>("alpaka")}) {
      std::vector<edm::EDGetTokenT<edm::ValueMap<float>>> v;
      for (auto const& n : kF)
        v.push_back(consumes<edm::ValueMap<float>>(edm::InputTag(label, n)));
      fTok_.push_back(v);
      posTok_.push_back(consumes<edm::ValueMap<GlobalPoint>>(edm::InputTag(label, "generalTrackmtdpos")));
      assocTok_.push_back(consumes<edm::ValueMap<int>>(edm::InputTag(label, "generalTrackassoc")));
    }
    std::fprintf(out_, "event,itrack,pt,eta,phi,nvalid,t_qual,t_pdg,t_vt,t_bx,t_evt");
    for (char const* p : {"L_", "N_"}) {
      for (auto const& s : kShort)
        std::fprintf(out_, ",%s%s", p, s.c_str());
      std::fprintf(out_, ",%sposx,%sposy,%sposz,%sassoc", p, p, p, p);
    }
    std::fprintf(out_, "\n");
  }
  ~MtdTimingCompare() override { std::fclose(out_); }

  void analyze(edm::Event const& ev, edm::EventSetup const&) override {
    auto const tracksH = ev.getHandle(tracksToken_);
    auto const r2sH = ev.getHandle(r2sToken_);
    std::vector<std::vector<edm::ValueMap<float> const*>> maps(2);
    std::vector<edm::ValueMap<GlobalPoint> const*> pos(2);
    std::vector<edm::ValueMap<int> const*> assoc(2);
    for (int k = 0; k < 2; ++k) {
      for (auto const& t : fTok_[k])
        maps[k].push_back(&ev.get(t));
      pos[k] = &ev.get(posTok_[k]);
      assoc[k] = &ev.get(assocTok_[k]);
    }
    for (size_t i = 0; i < tracksH->size(); ++i) {
      reco::TrackRef ref(tracksH, i);
      auto const& trk = (*tracksH)[i];
      double qual = -1, vt = 0;
      int pdg = 0, bx = -99, evt = -99;
      edm::RefToBase<reco::Track> rb(ref);
      if (r2sH.isValid()) {
        auto it = r2sH->find(rb);
        if (it != r2sH->end()) {
          for (auto const& tq : it->val)
            if (tq.second > qual) {
              qual = tq.second;
              auto const& tp = *tq.first;
              pdg = tp.g4Tracks().empty() ? 0 : tp.g4Tracks()[0].type();
              auto const& v = *tp.parentVertex();
              vt = v.position().t();
              bx = v.eventId().bunchCrossing();
              evt = v.eventId().event();
            }
        }
      }
      std::fprintf(out_,
                   "%llu,%zu,%.6g,%.6g,%.6g,%d,%.3f,%d,%.9g,%d,%d",
                   (unsigned long long)ev.id().event(),
                   i,
                   trk.pt(),
                   trk.eta(),
                   trk.phi(),
                   trk.numberOfValidHits(),
                   qual,
                   pdg,
                   vt,
                   bx,
                   evt);
      for (int k = 0; k < 2; ++k) {
        for (auto const* m : maps[k])
          std::fprintf(out_, ",%.9g", (*m)[ref]);
        auto const& p = (*pos[k])[ref];
        std::fprintf(out_, ",%.6g,%.6g,%.6g,%d", p.x(), p.y(), p.z(), (*assoc[k])[ref]);
      }
      std::fprintf(out_, "\n");
    }
  }

private:
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> r2sToken_;
  std::vector<std::vector<edm::EDGetTokenT<edm::ValueMap<float>>>> fTok_;
  std::vector<edm::EDGetTokenT<edm::ValueMap<GlobalPoint>>> posTok_;
  std::vector<edm::EDGetTokenT<edm::ValueMap<int>>> assocTok_;
  FILE* out_;
};

DEFINE_FWK_MODULE(MtdTimingCompare);
