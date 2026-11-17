#ifndef RecoVertex_PrimaryVertexProducer_interface_TrackForVertexFiller_h
#define RecoVertex_PrimaryVertexProducer_interface_TrackForVertexFiller_h

// Host-side filling of the TrackForVertexSoA used by the Alpaka DA in blocks: select the tracks as the
// standard producer does, compute their state at the beam line, sort them by z. Shared by
// PortableTrackSoAProducer and by the single packed input of the MTD + vertex chain.

#include <algorithm>
#include <cmath>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexSoA/interface/TrackForVertexSoA.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

namespace pvalpaka {

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

    static edm::ParameterSetDescription description() {
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
      return psd0;
    }

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

  // Fill one row if the track passes the selection; returns its weight, or -1 if it does not pass.
  template <typename TRow>
  double convertTrack(TRow out,
                      reco::TransientTrack const& in,
                      reco::BeamSpot const& bs,
                      TrackFilterParametersForVertexing const& fParams,
                      int32_t idx,
                      int32_t order) {
    double weight = -1;
    // stateAtBeamLine() returns by value, and so do the quantities derived from it: evaluate each of them once
    const auto stateAtBeamLine = in.stateAtBeamLine();
    const auto transverseIP = stateAtBeamLine.transverseImpactParameter();
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
      out.x() = position.x();
      out.y() = position.y();
      out.z() = position.z();
      out.px() = momentum.x();
      out.py() = momentum.y();
      out.pz() = momentum.z();
      out.weight() = weight;
      // the index in the reco::Track collection, to go back to it when converting the vertices
      out.ttIndex() = idx;
      out.dz2() = std::pow(in.track().dzError(), 2);
      // dz^2 + (bs*pt)^2*pz^2/pt^2 + vertexSize^2: accounts for the correlations and the vertex size
      double oneoverdz2 = (out.dz2()) +
                          ((std::pow(bs.BeamWidthX() * out.px(), 2)) + (std::pow(bs.BeamWidthY() * out.py(), 2))) *
                              std::pow(out.pz(), 2) / std::pow(momentum.perp2(), 2) +
                          std::pow(fParams.vertexSize, 2);
      oneoverdz2 = 1. / oneoverdz2;
      out.oneoverdz2() = oneoverdz2;
      out.dxy2AtIP() = std::pow(in.track().dxyError(), 2);
      out.dxy2() = std::pow(transverseIP.error(), 2);
      out.order() = order;
      out.isGood() = true;
    }
    return weight;
  }

  // Fill the view (capacity >= number of tracks) with the selected tracks sorted by z at the beam line.
  // Sets nT() and totweight(). Returns the number of selected tracks.
  inline int32_t fillTracksForVertexing(reco::TrackForVertexSoA::View tview,
                                        edm::Handle<reco::TrackCollection> const& tracks,
                                        reco::BeamSpot const& beamSpot,
                                        TransientTrackBuilder const& theB,
                                        TrackFilterParametersForVertexing const& fParams) {
    int32_t const tsize = tracks->size();
    tview.totweight() = 0;
    tview.nT() = 0;
    std::vector<reco::TransientTrack> t_tks = theB.build(tracks, beamSpot);
    // Sort by z at the beam line. stateAtBeamLine() returns by value, so evaluate z once per track.
    using ZType = std::decay_t<decltype(t_tks.front().stateAtBeamLine().trackStateAtPCA().position().z())>;
    std::vector<std::pair<ZType, int32_t>> sortedZAndIndex;
    sortedZAndIndex.reserve(tsize);
    for (int32_t idx = 0; idx < tsize; idx++)
      sortedZAndIndex.emplace_back(t_tks[idx].stateAtBeamLine().trackStateAtPCA().position().z(), idx);
    std::sort(sortedZAndIndex.begin(), sortedZAndIndex.end(), [](auto const& a, auto const& b) {
      return a.first < b.first;
    });
    int32_t nTrueTracks = 0;
    for (int32_t idx = 0; idx < tsize; idx++) {
      const int32_t trackIdx = sortedZAndIndex[idx].second;
      double weight = convertTrack(tview[nTrueTracks], t_tks[trackIdx], beamSpot, fParams, trackIdx, nTrueTracks);
      if (weight > 0) {
        nTrueTracks += 1;
        tview.nT() += 1;
        tview.totweight() += weight;
      }
    }
    return nTrueTracks;
  }

}  // namespace pvalpaka

#endif
