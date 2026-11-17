// Packs everything the device chain (MTD timing, then vertices) needs from one event into one host
// product (MtdEventInputHost): the tracks at their outermost hit, the MTD hits sorted into search
// cells, per-layer numbers measured on the hits, and the tracks selected for vertexing sorted by z.
// Device modules that consume it get it with a single host-to-device copy, made by the framework
// (no copy at all on the CPU backend). Only host work is done here.

#include <algorithm>
#include <cmath>
#include <vector>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/MTDDetId.h"
#include "DataFormats/MtdTimingSoA/interface/MtdEventInputHost.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CommonTopologies/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"
#include "RecoMTD/TrackExtenderAlpaka/interface/MtdBtlGeometryHost.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackForVertexFiller.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace ::mtdtiming;

  class MtdEventInputProducer : public global::EDProducer<> {
  public:
    explicit MtdEventInputProducer(edm::ParameterSet const& cfg)
        : EDProducer<>(cfg),
          tracksToken_{consumes(cfg.getParameter<edm::InputTag>("tracksSrc"))},
          hitsToken_{consumes(cfg.getParameter<edm::InputTag>("hitsSrc"))},
          beamSpotToken_{consumes(cfg.getParameter<edm::InputTag>("beamSpotSrc"))},
          geomToken_{esConsumes()},
          fieldToken_{esConsumes()},
          btlToken_{esConsumes()},
          ttbToken_{esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))},
          filter_{cfg.getParameter<edm::ParameterSet>("TkFilterParameters")},
          putToken_{produces()} {}

    void produce(edm::StreamID, device::Event& ev, device::EventSetup const& es) const override {
      auto const tracksHandle = ev.getHandle(tracksToken_);
      auto const& tracks = *tracksHandle;
      auto const& hits = ev.get(hitsToken_);
      auto const& geom = es.getData(geomToken_);
      auto const& field = es.getData(fieldToken_);
      auto const btl = es.getData(btlToken_).const_view().lookup();
      int32_t const nTracks = tracks.size();

      // ---------------------------------------------------------- hits: layer and search cell
      struct HitRef {
        int cell;
        int layer;
        GeomDet const* det;
        MTDTrackingRecHit const* hit;
      };
      std::vector<HitRef> refs;
      refs.reserve(hits.dataSize());
      float layerRef[kNLayers], layerSpread[kNLayers], maxHitErr[kNLayers];
      double sumRef[kNLayers] = {0, 0, 0, 0, 0};
      int nRef[kNLayers] = {0, 0, 0, 0, 0};
      for (int l = 0; l < kNLayers; ++l) {
        layerSpread[l] = 0.f;
        maxHitErr[l] = 0.f;
      }
      layerRef[0] = btl.radius();
      uint32_t const* lookupBegin = btl.rawId().data();
      auto const* lookupEnd = lookupBegin + btl.metadata().size();
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
          // for BTL the search cell is the module, as the legacy navigation sees it
          auto const it = std::lower_bound(lookupBegin, lookupEnd, id.rawId());
          if (it == lookupEnd || *it != id.rawId())
            continue;  // not part of the legacy navigation: the legacy can never find these hits
          btlCell = btl[it - lookupBegin].module();
        }
        for (auto const& hit : detSet) {
          auto const g = det->surface().toGlobal(hit.localPosition());
          refs.push_back({layer == 0 ? btlCell : CellGrid::cellOf(layer, g.x(), g.y(), g.z()), layer, det, &hit});
          auto const e = hit.localPositionError();
          maxHitErr[layer] = std::max(maxHitErr[layer], std::sqrt(std::max(e.xx(), e.yy())));
        }
      }
      for (int l = 1; l < kNLayers; ++l)
        layerRef[l] = nRef[l] > 0 ? sumRef[l] / nRef[l] : (l < 3 ? -300.f : 300.f);
      for (auto const& r : refs) {
        auto const& sp = r.det->surface().position();
        float const d = r.layer == 0 ? std::hypot(sp.x(), sp.y()) - layerRef[0] : sp.z() - layerRef[r.layer];
        layerSpread[r.layer] = std::max(layerSpread[r.layer], std::abs(d));
      }
      // stable sort by cell keeps the input order inside a cell: deterministic
      std::stable_sort(refs.begin(), refs.end(), [](HitRef const& a, HitRef const& b) { return a.cell < b.cell; });
      int32_t const nHits = refs.size();

      // ---------------------------------------------------------- one buffer for everything
      MtdEventInputHost out(ev.queue(), nTracks, nHits, CellGrid::kNCells + 1, kNLayers, nTracks);
      auto view = out.view();

      auto const& beamSpot = ev.get(beamSpotToken_);
      auto tv = view.tracks();
      for (int32_t i = 0; i < nTracks; ++i) {
        auto const& trk = tracks[i];
        auto row = tv[i];
        row.pt() = trk.pt();
        row.eta() = trk.eta();
        row.phi() = trk.phi();
        row.dzBS() = trk.dz(beamSpot.position());
        row.dxyBS() = trk.dxy(beamSpot.position());
        row.chi2() = trk.chi2();
        row.ndof() = trk.ndof();
        row.npixBarrel() = trk.hitPattern().numberOfValidPixelBarrelHits();
        row.npixEndcap() = trk.hitPattern().numberOfValidPixelEndcapHits();
        row.vx() = trk.vx();
        row.vy() = trk.vy();
        row.vz() = trk.vz();
        row.px() = trk.px();
        row.py() = trk.py();
        row.pz() = trk.pz();
        row.dzError() = trk.dzError();
        row.outerRadius() = row.outerZ() = 0.f;
        row.valid() = trk.extra().isAvailable() ? 1 : 0;
        row.charge() = trk.charge();
        row.x0() = trk.vx();
        row.y0() = trk.vy();
        row.z0() = trk.vz();
        row.ptPca() = trk.pt();
        row.pPca() = trk.p();
        row.cosThetaPca() = trk.pz() / trk.p();
        row.sigmaPPca() =
            std::sqrt(trk.covariance(reco::TrackBase::i_qoverp, reco::TrackBase::i_qoverp)) * trk.p() * trk.p();
        if (!row.valid()) {
          row.bz() = 0.f;
          row.outX() = row.outY() = row.outZ() = 0.f;
          row.outPx() = row.outPy() = row.outPz() = 0.f;
          row.outCov() = Vector15f::Zero();
          continue;
        }
        row.outerRadius() = trk.outerRadius();
        row.outerZ() = trk.outerZ();
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

      auto hv = view.hits();
      auto cv = view.cells();
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

      auto lv = view.layers();
      for (int l = 0; l < kNLayers; ++l) {
        lv[l].ref() = layerRef[l];
        lv[l].spread() = layerSpread[l];
        lv[l].maxHitErr() = maxHitErr[l];
      }

      // ---------------------------------------------------------- tracks for vertexing
      pvalpaka::fillTracksForVertexing(view.vtxTracks(), tracksHandle, beamSpot, es.getData(ttbToken_), filter_);

      ev.emplace(putToken_, std::move(out));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("tracksSrc", edm::InputTag("generalTracks"));
      desc.add<edm::InputTag>("hitsSrc", edm::InputTag("mtdTrackingRecHits"));
      desc.add<edm::InputTag>("beamSpotSrc", edm::InputTag("offlineBeamSpot"));
      desc.add<edm::ParameterSetDescription>("TkFilterParameters",
                                             pvalpaka::TrackFilterParametersForVertexing::description());
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    edm::EDGetTokenT<reco::TrackCollection> const tracksToken_;
    edm::EDGetTokenT<MTDTrackingDetSetVector> const hitsToken_;
    edm::EDGetTokenT<reco::BeamSpot> const beamSpotToken_;
    edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> const geomToken_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> const fieldToken_;
    edm::ESGetToken<MtdBtlGeometryHost, MTDRecoGeometryRecord> const btlToken_;
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> const ttbToken_;
    pvalpaka::TrackFilterParametersForVertexing const filter_;
    edm::EDPutTokenT<MtdEventInputHost> const putToken_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(MtdEventInputProducer);
