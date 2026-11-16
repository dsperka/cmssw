// MTD track timing, Alpaka version (no refit).
//
// Host part (this file): pack the tracks and the MTD hits into flat tables, sort the hits into
// cells, copy to the device. Device part (MtdTimingAlgo.dev.cc): find the best hit per layer and
// compute the timing. Output: one MtdTiming row per input track.

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/MTDDetId.h"
#include "DataFormats/MtdTimingSoA/interface/alpaka/MtdTimingDevice.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CommonTopologies/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"
#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"
#include "TrackingTools/DetLayers/interface/DetRod.h"
#include "TrackingTools/DetLayers/interface/PhiBorderFinder.h"
#include <unordered_map>
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "MtdTimingAlgo.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace ::mtdtiming;
  using ::mtdtiming::CellGrid;
  using ::mtdtiming::MatchParams;

  class MtdTimingProducer : public global::EDProducer<> {
  public:
    explicit MtdTimingProducer(edm::ParameterSet const& cfg)
        : EDProducer<>(cfg),
          tracksToken_{consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("tracksSrc"))},
          hitsToken_{consumes<MTDTrackingDetSetVector>(cfg.getParameter<edm::InputTag>("hitsSrc"))},
          geomToken_{esConsumes()},
          fieldToken_{esConsumes()},
          layersToken_{esConsumes()},
          outToken_{produces()} {
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
      auto const& tracks = ev.get(tracksToken_);
      auto const& hits = ev.get(hitsToken_);
      auto const& geom = es.getData(geomToken_);
      auto const& field = es.getData(fieldToken_);
      auto& queue = ev.queue();
      int32_t const nTracks = tracks.size();
      auto const tStart = std::chrono::steady_clock::now();

      // ---------------------------------------------------------- tracks
      PortableHostCollection<MtdTrackInSoA> tHost(queue, nTracks);
      auto tv = tHost.view();
      for (int32_t i = 0; i < nTracks; ++i) {
        auto const& trk = tracks[i];
        auto row = tv[i];
        row.valid() = trk.extra().isAvailable() ? 1 : 0;
        row.charge() = trk.charge();
        row.x0() = trk.vx();
        row.y0() = trk.vy();
        row.z0() = trk.vz();
        row.ptPca() = trk.pt();
        row.pPca() = trk.p();
        row.cosThetaPca() = trk.pz() / trk.p();
        row.sigmaPPca() = std::sqrt(trk.covariance(reco::TrackBase::i_qoverp, reco::TrackBase::i_qoverp)) * trk.p() * trk.p();
        if (!row.valid()) {
          row.bz() = 0.f;
          row.outX() = row.outY() = row.outZ() = 0.f;
          row.outPx() = row.outPy() = row.outPz() = 0.f;
          row.outCov() = Vector15f::Zero();
          continue;
        }
        auto const& op = trk.outerPosition();
        auto const& om = trk.outerMomentum();
        row.outX() = op.x();
        row.outY() = op.y();
        row.outZ() = op.z();
        row.outPx() = om.x();
        row.outPy() = om.y();
        row.outPz() = om.z();
        row.bz() = field.inTesla(GlobalPoint(op.x(), op.y(), op.z())).z();
        auto const& c = trk.extra()->outerStateCovariance();
        Vector15f cv;
        int k = 0;
        for (int a = 0; a < 5; ++a)
          for (int b = 0; b <= a; ++b)
            cv(k++) = c(a, b);
        row.outCov() = cv;
      }

      // ---------------------------------------------------------- BTL geometry seen by the legacy navigation
      // rods = rows of modules along z, sorted in phi (MTDTrayBarrelLayer); modules in rod order.
      auto const& layers = es.getData(layersToken_);
      auto const* btlLayer = static_cast<BarrelDetLayer const*>(layers.allBTLLayers().front());
      auto const& rods = btlLayer->components();
      std::unordered_map<uint32_t, int32_t> btlModuleIndex;
      int32_t nBtlModules = 0;
      for (auto const* rod : rods)
        nBtlModules += rod->basicComponents().size();
      PortableHostCollection<MtdPlaneSoA> rodHost(queue, int32_t(rods.size()));
      PortableHostCollection<MtdPlaneSoA> modHost(queue, nBtlModules);
      auto fillPlane = [](auto& view, int32_t i, Surface const& surf) {
        auto row = view[i];
        auto const& p = surf.position();
        auto const& r = surf.rotation();
        row.px() = p.x();
        row.py() = p.y();
        row.pz() = p.z();
        row.axX() = r.xx();
        row.axY() = r.xy();
        row.axZ() = r.xz();
        row.ayX() = r.yx();
        row.ayY() = r.yy();
        row.ayZ() = r.yz();
        row.nX() = r.zx();
        row.nY() = r.zy();
        row.nZ() = r.zz();
        row.halfX() = 0.5f * surf.bounds().width();
        row.halfY() = 0.5f * surf.bounds().length();
      };
      {
        int32_t m = 0;
        for (int32_t r = 0; r < int32_t(rods.size()); ++r) {
          auto rodView = rodHost.view();
          fillPlane(rodView, r, rods[r]->surface());
          auto rv = rodView[r];
          auto const dets = rods[r]->basicComponents();
          rv.first() = m;
          rv.count() = dets.size();
          for (auto const* det : dets) {
            auto modView = modHost.view();
            fillPlane(modView, m, det->surface());
            auto mv = modView[m];
            mv.first() = r;
            mv.count() = 0;
            mv.border() = 0.f;
            btlModuleIndex[det->geographicalId().rawId()] = m;
            ++m;
          }
        }
      }
      MatchParams params = params_;
      {
        auto const& cyl = btlLayer->specificSurface();
        params.layerRef[0] = cyl.radius();
        params.btlHalfLength = 0.5f * cyl.bounds().length();
        params.btlHalfThickness = 0.5f * cyl.bounds().thickness();
        params.btlNRods = rods.size();
        params.btlRodPhiStep = 2. * M_PI / rods.size();
        params.btlRodPhiOffset = rods.front()->position().phi() - 0.5 * params.btlRodPhiStep;
        std::vector<DetRod const*> detRods;
        for (auto const* rod : rods)
          detRods.push_back(static_cast<DetRod const*>(rod));
        PhiBorderFinder const bf(detRods);
        params.btlRodsOverlap = bf.isPhiOverlapping() ? 1 : 0;
        params.btlRodsPeriodic = bf.isPhiPeriodic() ? 1 : 0;
        for (int32_t r = 0; r < int32_t(rods.size()); ++r)
          rodHost.view()[r].border() = bf.phiBorders()[r];
      }

      // ---------------------------------------------------------- hits: count per cell, then fill sorted
      struct HitRef {
        int cell;
        int layer;
        GeomDet const* det;
        MTDTrackingRecHit const* hit;
      };
      std::vector<HitRef> refs;
      refs.reserve(hits.dataSize());
      params.debug = std::getenv("MTDT_DEBUGCAND") ? int(ev.id().event()) : 0;
      params.skipCandidates = std::getenv("MTDT_NOCAND") ? 1 : 0;
      double sumRef[kNLayers] = {0, 0, 0, 0, 0};
      int nRef[kNLayers] = {0, 0, 0, 0, 0};
      for (int l = 0; l < kNLayers; ++l) {
        params.layerSpread[l] = 0.f;
        params.maxHitErr[l] = 0.f;
      }
      for (auto const& detSet : hits) {
        DetId const id(detSet.detId());
        GeomDet const* det = geom.idToDet(id);
        if (det == nullptr)
          continue;
        MTDDetId const mid(id);
        int layer = 0;
        if (mid.mtdSubDetector() == MTDDetId::ETL) {
          ETLDetId const eid(id);
          layer = (eid.zside() > 0 ? 3 : 1) + (eid.nDisc() - 1);
        }
        auto const& sp = det->surface().position();
        sumRef[layer] += layer == 0 ? std::hypot(sp.x(), sp.y()) : sp.z();
        nRef[layer] += 1;
        int btlCell = -1;
        if (layer == 0) {
          auto const it = btlModuleIndex.find(id.rawId());
          if (it == btlModuleIndex.end())
            continue;  // not part of the legacy navigation: the legacy can never find these hits
          btlCell = it->second;
        }
        for (auto const& hit : detSet) {
          auto const g = det->surface().toGlobal(hit.localPosition());
          refs.push_back({layer == 0 ? btlCell : CellGrid::cellOf(layer, g.x(), g.y(), g.z()), layer, det, &hit});
          auto const e = hit.localPositionError();
          params.maxHitErr[layer] = std::max(params.maxHitErr[layer], std::sqrt(std::max(e.xx(), e.yy())));
        }
      }
      for (int l = 1; l < kNLayers; ++l)
        params.layerRef[l] = nRef[l] > 0 ? sumRef[l] / nRef[l] : (l < 3 ? -300.f : 300.f);
      for (auto const& r : refs) {
        auto const& sp = r.det->surface().position();
        float const d = r.layer == 0 ? std::hypot(sp.x(), sp.y()) - params.layerRef[0] : sp.z() - params.layerRef[r.layer];
        params.layerSpread[r.layer] = std::max(params.layerSpread[r.layer], std::abs(d));
      }
      // stable sort by cell keeps the input order inside a cell: deterministic
      std::stable_sort(refs.begin(), refs.end(), [](HitRef const& a, HitRef const& b) { return a.cell < b.cell; });


      int32_t const nHits = refs.size();
      PortableHostCollection<MtdHitInSoA> hHost(queue, nHits);
      PortableHostCollection<MtdCellSoA> cHost(queue, CellGrid::kNCells + 1);
      auto hv = hHost.view();
      auto cv = cHost.view();
      for (int32_t c = 0; c <= CellGrid::kNCells; ++c)
        cv[c].start() = 0;
      for (int32_t h = 0; h < nHits; ++h) {
        auto const& r = refs[h];
        auto const& surf = r.det->surface();
        auto const& rot = surf.rotation();
        auto const& sp = surf.position();
        auto const& hit = *r.hit;
        auto const lp = hit.localPosition();
        auto const le = hit.localPositionError();
        auto const g = surf.toGlobal(lp);
        auto row = hv[h];
        row.layer() = r.layer;
        row.detId() = hit.geographicalId().rawId();
        row.orig() = int32_t(r.hit - hits.data().data());
        row.planeX() = sp.x();
        row.planeY() = sp.y();
        row.planeZ() = sp.z();
        row.axX() = rot.xx();
        row.axY() = rot.xy();
        row.axZ() = rot.xz();
        row.ayX() = rot.yx();
        row.ayY() = rot.yy();
        row.ayZ() = rot.yz();
        row.nX() = rot.zx();
        row.nY() = rot.zy();
        row.nZ() = rot.zz();
        row.halfX() = 0.5f * surf.bounds().width();
        row.halfY() = 0.5f * surf.bounds().length();
        row.lx() = lp.x();
        row.ly() = lp.y();
        row.exx() = le.xx();
        row.exy() = le.xy();
        row.eyy() = le.yy();
        row.time() = hit.time();
        row.timeError() = hit.timeError();
        row.gx() = g.x();
        row.gy() = g.y();
        row.gz() = g.z();
        cv[r.cell + 1].start() += 1;
      }
      for (int32_t c = 0; c < CellGrid::kNCells; ++c)
        cv[c + 1].start() += cv[c].start();

      // ---------------------------------------------------------- to the device, run, output
      mtdtiming::TrackInDevice tDev(queue, nTracks);
      mtdtiming::HitInDevice hDev(queue, nHits);
      mtdtiming::CellDevice cDev(queue, CellGrid::kNCells + 1);
      alpaka::memcpy(queue, tDev.buffer(), tHost.buffer());
      alpaka::memcpy(queue, hDev.buffer(), hHost.buffer());
      alpaka::memcpy(queue, cDev.buffer(), cHost.buffer());
      mtdtiming::PlaneDevice rodDev(queue, int32_t(rods.size()));
      mtdtiming::PlaneDevice modDev(queue, nBtlModules);
      alpaka::memcpy(queue, rodDev.buffer(), rodHost.buffer());
      alpaka::memcpy(queue, modDev.buffer(), modHost.buffer());
      mtdtiming::RecordDevice records(queue, nTracks * kNLayers);
      mtdtiming::MtdTimingDevice out(queue, nTracks);
      auto const tPacked = std::chrono::steady_clock::now();
      mtdtiming::runMtdTiming(queue, params, tDev, hDev, cDev, rodDev, modDev, records, out);
      if (std::getenv("MTDT_TIMING")) {
        alpaka::wait(queue);
        auto const tDone = std::chrono::steady_clock::now();
        printf("TIMING event %llu tracks %d hits %d pack %.1f ms kernels %.1f ms\n", (unsigned long long)ev.id().event(), nTracks, nHits,
               std::chrono::duration<double, std::milli>(tPacked - tStart).count(),
               std::chrono::duration<double, std::milli>(tDone - tPacked).count());
      }
      ev.emplace(outToken_, std::move(out));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("tracksSrc", edm::InputTag("generalTracks"));
      desc.add<edm::InputTag>("hitsSrc", edm::InputTag("mtdTrackingRecHits"));
      edm::ParameterSetDescription m;  // same names and defaults as the legacy MTDHitMatcher PSet
      m.add<double>("estimatorMaxChi2", 500.);
      m.add<double>("estimatorMaxNSigma", 10.);
      m.add<double>("btlChi2Cut", 50.);
      m.add<double>("btlTimeChi2Cut", 10.);
      m.add<double>("etlChi2Cut", 50.);
      m.add<double>("etlTimeChi2Cut", 10.);
      m.add<double>("bsTimeSpread", 0.2);
      desc.add<edm::ParameterSetDescription>("MTDHitMatcher", m);
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    edm::EDGetTokenT<reco::TrackCollection> const tracksToken_;
    edm::EDGetTokenT<MTDTrackingDetSetVector> const hitsToken_;
    edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> const geomToken_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> const fieldToken_;
    edm::ESGetToken<MTDDetLayerGeometry, MTDRecoGeometryRecord> const layersToken_;
    device::EDPutToken<mtdtiming::MtdTimingDevice> const outToken_;
    MatchParams params_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(MtdTimingProducer);
