#include <cmath>

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "vdt/vdtMath.h"

#include "VertexTimeAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using namespace cms::alpakatools;

  namespace {
    // the host algorithm uses vdt::fast_exp; use the same on the CPU so that the results are identical
    ALPAKA_FN_HOST_ACC inline double fastExp(double x) {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      return std::exp(x);
#else
      return vdt::fast_exp(x);
#endif
    }
  }  // namespace

  // one thread per vertex
  class VertexTimeKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  VertexTimeParameters params,
                                  ::reco::TrackForVertexSoA::ConstView tracksInBlocks,
                                  ::reco::VertexTrackSoA::ConstView vertexTracks,
                                  ::mtdtiming::MtdTimingSoA::ConstView timing,
                                  ::reco::VertexSoA::View vertices) const {
      int const nV = vertices[0].nV();
      double const a[3] = {params.probPion, params.probKaon, params.probProton};
      for (int i : uniform_elements(acc, nV)) {
        int const iV = vertices[i].order();
        vertices[iV].t() = 0.f;
        vertices[iV].errt() = -1.f;
        if (not vertices[iV].isGood())
          continue;
        int const first = vertices[iV].first_track();
        int const nTrk = vertices[iV].ntracks();

        // the time hypotheses of a track of the vertex, -1 weight if the track is not used
        auto trackInfo = [&](int iT, double& w, double tHyp[3], double eHyp[3]) {
          w = vertexTracks[first + iT].track_weight();
          if (!(w > params.minTrackVtxWeight)) {
            w = -1;
            return;
          }
          auto const r = timing[tracksInBlocks[vertexTracks[first + iT].track_id()].ttIndex()];
          if (params.useMVAVtxTime && !(r.mva() >= params.minTrackTimeQuality)) {
            w = -1;
            return;
          }
          float const e = r.sigmatmtd();
          eHyp[0] = std::sqrt(e * e + r.sigmaTofPi() * r.sigmaTofPi());
          eHyp[1] = std::sqrt(e * e + r.sigmaTofK() * r.sigmaTofK());
          eHyp[2] = std::sqrt(e * e + r.sigmaTofP() * r.sigmaTofP());
          tHyp[0] = r.tmtd() - r.tofPi();
          tHyp[1] = r.tmtd() - r.tofK();
          tHyp[2] = r.tmtd() - r.tofP();
        };

        double tsum = 0, wsum = 0, w2sum = 0;
        for (int iT = 0; iT < nTrk; ++iT) {
          double w, tHyp[3], eHyp[3];
          trackInfo(iT, w, tHyp, eHyp);
          if (w < 0)
            continue;
          for (int j = 0; j < 3; ++j) {
            double const wgt = w / (eHyp[j] * eHyp[j]);
            wsum += wgt * a[j];
            tsum += wgt * a[j] * tHyp[j];
          }
        }
        if (!(wsum > 0))
          continue;

        double t0 = tsum / wsum;
        double beta = 1. / params.Tstart;
        int nit = 0;
        constexpr int maxIterations = 100;
        while ((nit++) < maxIterations) {
          tsum = 0;
          wsum = 0;
          w2sum = 0;
          for (int iT = 0; iT < nTrk; ++iT) {
            double w, tHyp[3], eHyp[3];
            trackInfo(iT, w, tHyp, eHyp);
            if (w < 0)
              continue;
            double e[3] = {0, 0, 0};
            constexpr double cutOff = 4.5;
            double Z = fastExp(-beta * cutOff);  // outlier term exp(-beta * 0.5 * 3 * 3)
            for (int j = 0; j < 3; j++) {
              double const tpull = (tHyp[j] - t0) / eHyp[j];
              e[j] = fastExp(-0.5 * beta * tpull * tpull);
              Z += a[j] * e[j];
            }
            double wsumTrk = 0, wsumSigmaTrk = 0;
            for (int j = 0; j < 3; j++) {
              double const wt = a[j] * e[j] / Z;
              double const ww = wt * w / (eHyp[j] * eHyp[j]);
              wsumTrk += ww;
              wsumSigmaTrk += ww * eHyp[j];
              tsum += ww * tHyp[j];
            }
            wsum += wsumTrk;
            w2sum += wsumSigmaTrk * wsumSigmaTrk;
          }
          if (wsum < 1e-10)
            break;  // failed
          float const vtxTime = tsum / wsum;
          if ((std::abs(vtxTime - t0) < 1e-4 / std::sqrt(beta)) and beta >= 1.) {
            vertices[iV].t() = vtxTime;
            vertices[iV].errt() = float(std::sqrt(w2sum) / wsum);
            break;
          }
          if ((std::abs(vtxTime - t0) < 1e-3) and beta < 1.)
            beta = std::min(1., beta / params.coolingFactor);
          t0 = vtxTime;
        }
      }
    }
  };

  void vertexTimes(Queue& queue,
                   VertexTimeParameters const& params,
                   ::reco::TrackForVertexSoA::ConstView tracksInBlocks,
                   ::reco::VertexTrackSoA::ConstView vertexTracks,
                   ::mtdtiming::MtdTimingSoA::ConstView timing,
                   ::reco::VertexSoA::View vertices) {
    constexpr int threads = 64;
    alpaka::exec<Acc1D>(queue,
                        make_workdiv<Acc1D>(divide_up_by(1024, threads), threads),
                        VertexTimeKernel{},
                        params,
                        tracksInBlocks,
                        vertexTracks,
                        timing,
                        vertices);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
