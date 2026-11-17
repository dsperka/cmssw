// TOF particle identification on the device: the same algorithm as RecoMTD/TimingIDTools TOFPIDProducer,
// one thread per track, looping over the vertices.

#include <cmath>
#include <limits>

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "TofPidAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming {
  using namespace cms::alpakatools;

  namespace {
    // what TOFPIDProducer reads from a reco::Vertex made by RecoVertexFromSoAProducer
    struct VertexInfo {
      double x, y, z, t, tError;
    };

    ALPAKA_FN_ACC inline VertexInfo vertexInfo(::reco::VertexSoA::ConstView vertices, int iV) {
      auto const v = vertices[iV];
      // reco::Vertex keeps the time as a double and its covariance as floats: tError = sqrt(float(errt^2))
      double const tError = v.errt() >= 0 ? std::sqrt(double(float(double(v.errt()) * double(v.errt())))) : 0.;
      double const t = v.errt() >= 0 ? double(v.t()) : 0.;
      return VertexInfo{v.x(), v.y(), v.z(), t, tError};
    }

    // reco::TrackBase::dz(point)
    ALPAKA_FN_ACC inline double dzTo(::mtdtiming::MtdTrackInSoA::ConstView::const_element t, VertexInfo const& v) {
      double const pt = std::sqrt(t.px() * t.px() + t.py() * t.py());
      return (t.vz() - v.z) - ((t.vx() - v.x) * t.px() + (t.vy() - v.y) * t.py()) / pt * t.pz() / pt;
    }
  }  // namespace

  // the vertex in which each track has a weight above 0.5 (-1 if none)
  class TrackVertexKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  ::reco::VertexSoA::ConstView vertices,
                                  ::reco::VertexTrackSoA::ConstView vertexTracks,
                                  ::reco::TrackForVertexSoA::ConstView tracksInBlocks,
                                  int32_t* trackVertex,
                                  int32_t nTracks) const {
      for (int32_t i : uniform_elements(acc, nTracks))
        trackVertex[i] = -1;
    }
  };

  class AssignTracksKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  ::reco::VertexSoA::ConstView vertices,
                                  ::reco::VertexTrackSoA::ConstView vertexTracks,
                                  ::reco::TrackForVertexSoA::ConstView tracksInBlocks,
                                  int32_t* trackVertex) const {
      int const nV = vertices[0].nV();
      for (int i : uniform_elements(acc, nV)) {
        int const iV = vertices[i].order();
        if (not vertices[iV].isGood())
          continue;
        int const first = vertices[iV].first_track();
        for (int iT = 0; iT < vertices[iV].ntracks(); ++iT) {
          if (vertexTracks[first + iT].track_weight() > 0.5f)
            alpaka::atomicMax(
                acc, &trackVertex[tracksInBlocks[vertexTracks[first + iT].track_id()].ttIndex()], iV, alpaka::hierarchy::Grids{});
        }
      }
    }
  };

  class TofPidKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  TofPidParameters p,
                                  ::mtdtiming::MtdTrackInSoA::ConstView tracks,
                                  ::mtdtiming::MtdTimingSoA::ConstView timing,
                                  ::reco::VertexSoA::ConstView vertices,
                                  int32_t const* trackVertex,
                                  ::mtdtiming::TofPidSoA::View out) const {
      int const nV = vertices[0].nV();
      for (int32_t itrack : uniform_elements(acc, tracks.metadata().size())) {
        auto const track = tracks[itrack];
        auto const r = timing[itrack];
        float t0 = r.t0();
        float t0safe = t0;
        float sigmat0safe = r.sigmat0();
        float sigmatmtd = (r.sigmatmtd() > 0. && p.fixedT0Error > 0.) ? float(p.fixedT0Error) : r.sigmatmtd();
        float sigmat0 = sigmatmtd;
        float const sigmatofpi = r.sigmaTofPi();
        float const sigmatofk = r.sigmaTofK();
        float const sigmatofp = r.sigmaTofP();
        float prob_pi = -1.;
        float prob_k = -1.;
        float prob_p = -1.;
        float const trackMVAQual = r.mva();
        if (sigmat0 > 0. && (!p.MVASel || (p.MVASel && trackMVAQual >= p.minTrackTimeQuality))) {
          double const rsigmazsq = 1. / track.dzError() / track.dzError();
          double const rsigmat[3] = {1. / std::sqrt(sigmatmtd * sigmatmtd + sigmatofpi * sigmatofpi),
                                     1. / std::sqrt(sigmatmtd * sigmatmtd + sigmatofk * sigmatofk),
                                     1. / std::sqrt(sigmatmtd * sigmatmtd + sigmatofp * sigmatofp)};
          // nominal vertex: the one the track belongs to, else the closest in z (or in z and t)
          int vtxidx = trackVertex[itrack];
          int vtxidxmindz = -1;
          int vtxidxminchisq = -1;
          double mindz = p.maxDz;
          double minchisq = std::numeric_limits<double>::max();
          if (vtxidx < 0) {
            for (int i = 0; i < nV; ++i) {
              int const iV = vertices[i].order();
              if (not vertices[iV].isGood())
                continue;
              auto const vtx = vertexInfo(vertices, iV);
              double const dz = std::abs(dzTo(track, vtx));
              if (dz < mindz) {
                mindz = dz;
                vtxidxmindz = iV;
              }
              if (vtx.tError > 0. && vtx.tError < p.vtxMaxSigmaT) {
                double const dt = std::abs(t0 - vtx.t);
                double const dtsig = dt * rsigmat[0];
                double const chisq = dz * dz * rsigmazsq + dtsig * dtsig;
                if (dz < p.maxDz && dtsig < p.maxDtSignificance && chisq < minchisq) {
                  minchisq = chisq;
                  vtxidxminchisq = iV;
                }
              }
            }
            auto timed = [&](int iV) {
              double const e = vertexInfo(vertices, iV).tError;
              return e > 0. && e < p.vtxMaxSigmaT;
            };
            if (vtxidxmindz >= 0 && !timed(vtxidxmindz))
              vtxidx = vtxidxmindz;
            else if (vtxidxminchisq >= 0)
              vtxidx = vtxidxminchisq;
            else if (vtxidxmindz >= 0)
              vtxidx = vtxidxmindz;
          }
          VertexInfo vtxnom{};
          if (vtxidx >= 0)
            vtxnom = vertexInfo(vertices, vtxidx);
          if (vtxidx >= 0 && vtxnom.tError > 0. && vtxnom.tError < p.vtxMaxSigmaT) {
            double const dznom = std::abs(dzTo(track, vtxnom));
            double const dtnom = std::abs(t0 - vtxnom.t);
            double const dtsignom = dtnom * rsigmat[0];
            double const chisqnom = dznom * dznom * rsigmazsq + dtsignom * dtsignom;
            double t0_best = t0;
            if (dtsignom < p.maxDtSignificance) {
              sigmat0safe = 1. / rsigmat[0];
              sigmat0 = sigmat0safe;
            }
            double const tmtd = r.tmtd();
            double const t0_k = tmtd - r.tofK();
            double const t0_p = tmtd - r.tofP();
            double chisqmin = chisqnom;
            double const chisqmin_pi = chisqnom;
            double chisqmin_k = std::numeric_limits<double>::max();
            double chisqmin_p = std::numeric_limits<double>::max();
            for (int i = 0; i < nV; ++i) {
              int const iV = vertices[i].order();
              if (not vertices[iV].isGood())
                continue;
              if (!p.vertexReassignment && iV != vtxidx)
                continue;
              auto const vtx = vertexInfo(vertices, iV);
              if (!(vtx.tError > 0. && vtx.tError < p.vtxMaxSigmaT))
                continue;
              double const dz = std::abs(dzTo(track, vtx));
              if (dz >= p.maxDz)
                continue;
              double const chisqdz = dz * dz * rsigmazsq;
              double const dt_k = std::abs(t0_k - vtx.t);
              double const dtsig_k = dt_k * rsigmat[1];
              double const chisq_k = chisqdz + dtsig_k * dtsig_k;
              if (dtsig_k < p.maxDtSignificance && chisq_k < chisqmin_k)
                chisqmin_k = chisq_k;
              double const dt_p = std::abs(t0_p - vtx.t);
              double const dtsig_p = dt_p * rsigmat[2];
              double const chisq_p = chisqdz + dtsig_p * dtsig_p;
              if (dtsig_p < p.maxDtSignificance && chisq_p < chisqmin_p)
                chisqmin_p = chisq_p;
              if (dtsig_k < p.maxDtSignificance && chisq_k < chisqmin) {
                chisqmin = chisq_k;
                t0_best = t0_k;
                t0safe = t0_k;
                sigmat0safe = 1. / rsigmat[1];
                sigmat0 = sigmat0safe;
              }
              if (dtsig_p < p.maxDtSignificance && chisq_p < chisqmin) {
                chisqmin = chisq_p;
                t0_best = t0_p;
                t0safe = t0_p;
                sigmat0safe = 1. / rsigmat[2];
                sigmat0 = sigmat0safe;
              }
            }
            double const rawprob_pi = p.probPion * std::exp(-0.5 * chisqmin_pi);
            double const rawprob_k = p.probKaon * std::exp(-0.5 * chisqmin_k);
            double const rawprob_p = p.probProton * std::exp(-0.5 * chisqmin_p);
            double const normprob = 1. / (rawprob_pi + rawprob_k + rawprob_p);
            prob_pi = rawprob_pi * normprob;
            prob_k = rawprob_k * normprob;
            prob_p = rawprob_p * normprob;
            double const prob_heavy = 1. - prob_pi;
            if (prob_heavy > p.minProbHeavy)
              t0 = t0_best;
          }
        }
        auto o = out[itrack];
        o.t0() = t0;
        o.sigmat0() = sigmat0;
        o.t0safe() = t0safe;
        o.sigmat0safe() = sigmat0safe;
        o.probPi() = prob_pi;
        o.probK() = prob_k;
        o.probP() = prob_p;
      }
    }
  };

  void runTofPid(Queue& queue,
                 TofPidParameters const& params,
                 ::mtdtiming::MtdTrackInSoA::ConstView tracks,
                 ::mtdtiming::MtdTimingSoA::ConstView timing,
                 ::reco::VertexSoA::ConstView vertices,
                 ::reco::VertexTrackSoA::ConstView vertexTracks,
                 ::reco::TrackForVertexSoA::ConstView tracksInBlocks,
                 TofPidDevice& out) {
    int32_t const nTracks = tracks.metadata().size();
    if (nTracks == 0)
      return;
    auto trackVertex = make_device_buffer<int32_t[]>(queue, nTracks);
    constexpr uint32_t threads = 64;
    auto const wdTracks = make_workdiv<Acc1D>(divide_up_by(nTracks, threads), threads);
    alpaka::exec<Acc1D>(
        queue, wdTracks, TrackVertexKernel{}, vertices, vertexTracks, tracksInBlocks, trackVertex.data(), nTracks);
    alpaka::exec<Acc1D>(queue,
                        make_workdiv<Acc1D>(divide_up_by(1024, threads), threads),
                        AssignTracksKernel{},
                        vertices,
                        vertexTracks,
                        tracksInBlocks,
                        trackVertex.data());
    alpaka::exec<Acc1D>(
        queue, wdTracks, TofPidKernel{}, params, tracks, timing, vertices, trackVertex.data(), out.view());
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming
