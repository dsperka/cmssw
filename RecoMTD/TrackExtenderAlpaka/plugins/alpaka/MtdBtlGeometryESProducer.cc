// BTL geometry for the Alpaka MTD timing, as the legacy navigation sees it (MTDTrayBarrelLayer):
// "rods" (rows of modules along z, sorted in phi) and their modules (in rod order), plus a table of
// the modules sorted by detector id. Made once per geometry; the framework copies it to the device
// once, not every event.

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <vector>

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ModuleFactory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"
#include "RecoMTD/TrackExtenderAlpaka/interface/MtdBtlGeometryHost.h"
#include "RecoMTD/TrackExtenderAlpaka/interface/alpaka/MtdBtlGeometryDevice.h"
#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"
#include "TrackingTools/DetLayers/interface/DetRod.h"
#include "TrackingTools/DetLayers/interface/PhiBorderFinder.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class MtdBtlGeometryESProducer : public ESProducer {
  public:
    explicit MtdBtlGeometryESProducer(edm::ParameterSet const& cfg) : ESProducer(cfg) {
      auto cc = setWhatProduced(this);
      layersToken_ = cc.consumes();
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      descriptions.addWithDefaultLabel(desc);
    }

    std::unique_ptr<::mtdtiming::MtdBtlGeometryHost> produce(MTDRecoGeometryRecord const& record) {
      auto const& layers = record.get(layersToken_);
      auto const* btlLayer = static_cast<BarrelDetLayer const*>(layers.allBTLLayers().front());
      auto const& rods = btlLayer->components();
      int32_t const nRods = rods.size();
      int32_t nModules = 0;
      for (auto const* rod : rods)
        nModules += rod->basicComponents().size();

      auto product = std::make_unique<::mtdtiming::MtdBtlGeometryHost>(cms::alpakatools::host(), nRods, nModules, nModules);
      auto rodView = product->view().rods();
      auto modView = product->view().modules();
      auto lookup = product->view().lookup();

      auto fillPlane = [](auto row, Surface const& surf) {
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

      std::vector<uint32_t> rawIds(nModules);
      int32_t m = 0;
      for (int32_t r = 0; r < nRods; ++r) {
        fillPlane(rodView[r], rods[r]->surface());
        auto const dets = rods[r]->basicComponents();
        rodView[r].first() = m;
        rodView[r].count() = dets.size();
        for (auto const* det : dets) {
          fillPlane(modView[m], det->surface());
          modView[m].first() = r;
          modView[m].count() = 0;
          modView[m].border() = 0.f;
          rawIds[m] = det->geographicalId().rawId();
          ++m;
        }
      }

      // modules sorted by detector id, for the hit -> module lookup on the host
      std::vector<int32_t> order(nModules);
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](int32_t a, int32_t b) { return rawIds[a] < rawIds[b]; });
      for (int32_t i = 0; i < nModules; ++i) {
        lookup[i].rawId() = rawIds[order[i]];
        lookup[i].module() = order[i];
      }

      // layer cylinder and rod bin finder (PeriodicBinFinderInPhi or the general one on the borders)
      auto const& cyl = btlLayer->specificSurface();
      lookup.radius() = cyl.radius();
      lookup.halfLength() = 0.5f * cyl.bounds().length();
      lookup.halfThickness() = 0.5f * cyl.bounds().thickness();
      lookup.rodPhiStep() = 2. * M_PI / nRods;
      lookup.rodPhiOffset() = rods.front()->position().phi() - 0.5 * lookup.rodPhiStep();
      std::vector<DetRod const*> detRods;
      for (auto const* rod : rods)
        detRods.push_back(static_cast<DetRod const*>(rod));
      PhiBorderFinder const bf(detRods);
      lookup.rodsOverlap() = bf.isPhiOverlapping() ? 1 : 0;
      lookup.rodsPeriodic() = bf.isPhiPeriodic() ? 1 : 0;
      for (int32_t r = 0; r < nRods; ++r)
        rodView[r].border() = bf.phiBorders()[r];
      return product;
    }

  private:
    edm::ESGetToken<MTDDetLayerGeometry, MTDRecoGeometryRecord> layersToken_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_EVENTSETUP_ALPAKA_MODULE(MtdBtlGeometryESProducer);
