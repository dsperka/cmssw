#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "RecoMTD/TrackExtenderAlpaka/interface/MtdTimingMath.h"

#include "MtdTimingAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming {
  using namespace cms::alpakatools;
  using namespace ::mtdtiming;

  namespace {
    ALPAKA_FN_ACC inline State outerState(MtdTrackInSoA::ConstView::const_element t) {
      double const px = t.outPx(), py = t.outPy(), pz = t.outPz();
      double const p = std::sqrt(px * px + py * py + pz * pz);
      return State{{t.outX(), t.outY(), t.outZ()}, {px / p, py / p, pz / p}, p, t.charge(), t.bz()};
    }

    ALPAKA_FN_ACC inline float trackerPathOf(MtdTrackInSoA::ConstView::const_element t) {
      return float(trackerPath(t.x0(),
                               t.y0(),
                               t.z0(),
                               t.ptPca(),
                               t.cosThetaPca(),
                               V3{t.outX(), t.outY(), t.outZ()},
                               V3{t.outPx(), t.outPy(), t.outPz()},
                               t.bz()));
    }

    ALPAKA_FN_ACC inline TofInput tofInputOf(MtdTrackInSoA::ConstView::const_element t, float gap) {
      float const p2Out = t.outPx() * t.outPx() + t.outPy() * t.outPy() + t.outPz() * t.outPz();
      return TofInput{trackerPathOf(t), gap, t.pPca() * t.pPca(), p2Out, t.sigmaPPca()};
    }
  }  // namespace

  // ------------------------------------------------------------------ kernel 1
  class MatchHitsKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  MtdTrackInSoA::ConstView tracks,
                                  MtdHitInSoA::ConstView hits,
                                  MtdCellSoA::ConstView cells,
                                  MtdPlaneSoA::ConstView rods,
                                  MtdPlaneSoA::ConstView mods,
                                  MtdRecordSoA::View records,
                                  MatchParams params) const {
      int32_t const nRecords = tracks.metadata().size() * kNLayers;
      for (int32_t r : uniform_elements(acc, nRecords)) {
        int32_t const itrk = r / kNLayers;
        int const layer = r % kNLayers;
        records[r].hit() = -1;
        records[r].est() = -1.f;
        records[r].timeChi2() = -1.f;
        records[r].path() = 0.f;

        auto const t = tracks[itrk];
        if (!t.valid())
          continue;
        // ETL: only the two discs on the same side as the outer state (legacy rule)
        if (layer > 0 && (layer >= 3) != (t.outZ() > 0.f))
          continue;

        State const st = outerState(t);
        float cov[15];
        for (int k = 0; k < 15; ++k)
          cov[k] = t.outCov()(k);

        // ---- candidate evaluation, common to BTL and ETL: keep the best hit by the legacy ordering
        TofInput tin = tofInputOf(t, 0.f);
        int32_t best = -1;
        float bestEst = 0.f, bestTchi2 = 0.f, bestPath = 0.f;
        auto evalHit = [&](int32_t h, double tx, double ty, double txx, double txy, double tyy, double s) {
          auto const hit = hits[h];
          float const est =
              float(chi2Local(hit.lx() - tx, hit.ly() - ty, txx + hit.exx(), txy + hit.exy(), tyy + hit.eyy()));
          if (!(est <= params.estimatorMaxChi2))
            return;
          tin.gapPath = float(std::abs(s));
          float const tchi2 = timeChi2(tin, hit.time(), hit.timeError(), params.bsTimeSpread);
          if (best < 0 || betterHit(est, tchi2, bestEst, bestTchi2)) {
            best = h;
            bestEst = est;
            bestTchi2 = tchi2;
            bestPath = float(s);
          }
        };
        double const nSig = params.nSigmaWindow;  // 10, the legacy estimator's nSigma

        if (layer == 0) {
          // ================= BTL: the legacy navigation (MTDTrayBarrelLayer + MTDDetTray) =================
          // 1) BarrelDetLayer::compatible: cross the layer cylinder; z inside the half-length, or inside
          //    half-length + half-thickness * |dz|/dperp + 3 sigma(z)
          double sc;
          if (!crossCylinder(st, params.layerRef[0], sc))
            continue;
          V3 cpos, cdir;
          moveAlong(st, sc, cpos, cdir);
          if (std::abs(cpos.z) >= params.btlHalfLength) {
            double const rr = std::sqrt(cpos.x * cpos.x + cpos.y * cpos.y);
            V3 const n{cpos.x / rr, cpos.y / rr, 0.};
            double lx, ly, cxx, cxy, cyy, ps;
            if (!localPositionAndError(st, cov, cpos, n, V3{-n.y, n.x, 0.}, V3{0., 0., 1.}, lx, ly, cxx, cxy, cyy, ps))
              continue;
            double const dperp = std::sqrt(cdir.x * cdir.x + cdir.y * cdir.y);
            double const dz = params.btlHalfLength + params.btlHalfThickness * std::abs(cdir.z) / dperp + 3. * std::sqrt(cyy);
            if (std::abs(cpos.z) >= dz)
              continue;
          }
          // a module is "compatible" if the crossing is inside its bounds widened by nSig track sigmas;
          // then all its hits are candidates. Returns the local x and errors of the crossing.
          auto visitModule = [&](int32_t m, double& outX, double& outSx, double& outSy) -> bool {
            auto const md = mods[m];
            V3 const P{md.px(), md.py(), md.pz()}, N{md.nX(), md.nY(), md.nZ()};
            V3 const AX{md.axX(), md.axY(), md.axZ()}, AY{md.ayX(), md.ayY(), md.ayZ()};
            double tx, ty, txx, txy, tyy, s;
            if (!localPositionAndError(st, cov, P, N, AX, AY, tx, ty, txx, txy, tyy, s))
              return false;
            bool const inside = (std::abs(tx) < md.halfX() && std::abs(ty) < md.halfY()) ||
                                (std::abs(tx) < md.halfX() + nSig * std::sqrt(txx) &&
                                 std::abs(ty) < md.halfY() + nSig * std::sqrt(tyy));
            if (!inside)
              return false;
            outX = tx;
            outSx = std::sqrt(txx);
            outSy = std::sqrt(tyy);
            for (int32_t h = cells[m].start(); h < cells[m + 1].start(); ++h)
              evalHit(h, tx, ty, txx, txy, tyy, s);
            return true;
          };
          // MTDDetTray::compatibleDets: rod plane within bounds (nSig), closest module in z, then the
          // neighbours on each side while within reach and compatible
          auto visitRod = [&](int32_t ir, double& firstX, double& firstSx) -> bool {
            auto const rd = rods[ir];
            V3 const P{rd.px(), rd.py(), rd.pz()}, N{rd.nX(), rd.nY(), rd.nZ()};
            V3 const AX{rd.axX(), rd.axY(), rd.axZ()}, AY{rd.ayX(), rd.ayY(), rd.ayZ()};
            double rx, ry, rxx, rxy, ryy, rs;
            if (!localPositionAndError(st, cov, P, N, AX, AY, rx, ry, rxx, rxy, ryy, rs))
              return false;
            bool const inside = (std::abs(rx) < rd.halfX() && std::abs(ry) < rd.halfY()) ||
                                (std::abs(rx) < rd.halfX() + nSig * std::sqrt(rxx) &&
                                 std::abs(ry) < rd.halfY() + nSig * std::sqrt(ryy));
            if (!inside)
              return false;
            V3 start, sdir;
            moveAlong(st, rs, start, sdir);
            int32_t const first = rd.first(), n = rd.count();
            // GenericBinFinderInZ: borders half way between consecutive module z positions
            int bin = 0;
            for (int i = 0; i < n - 1; ++i)
              if (start.z > 0.5 * (mods[first + i].pz() + mods[first + i + 1].pz()))
                bin = i + 1;
            double fx, fsx, fsy;
            if (!visitModule(first + bin, fx, fsx, fsy))
              return false;  // closest not compatible: the neighbours are not tried
            firstX = fx;
            firstSx = fsx;
            double const maxDy = nSig * fsy;
            double const halfLen = mods[first + bin].halfY();
            double dummy1, dummy2, dummy3;
            for (int i = bin + 1; i < n; ++i) {
              auto const md = mods[first + i];
              double const y = dot(sub(start, V3{md.px(), md.py(), md.pz()}), V3{md.ayX(), md.ayY(), md.ayZ()});
              if (!(std::abs(y) < halfLen + maxDy) || !visitModule(first + i, dummy1, dummy2, dummy3))
                break;
            }
            for (int i = bin - 1; i >= 0; --i) {
              auto const md = mods[first + i];
              double const y = dot(sub(start, V3{md.px(), md.py(), md.pz()}), V3{md.ayX(), md.ayY(), md.ayZ()});
              if (!(std::abs(y) < halfLen + maxDy) || !visitModule(first + i, dummy1, dummy2, dummy3))
                break;
            }
            return true;
          };
          // MTDTrayBarrelLayer::compatibleDets: the closest rod in phi; the next rod only if the track
          // is within 3 sigma of the closest rod's edge (or the closest rod gave nothing)
          int32_t const nRods = params.btlNRods;
          double const phi = std::atan2(cpos.y, cpos.x);
          int32_t closest = -1;
          if (params.btlRodsPeriodic) {  // PeriodicBinFinderInPhi
            double tmp = std::fmod(phi - params.btlRodPhiOffset, 2. * M_PI) / params.btlRodPhiStep;
            if (tmp < 0.)
              tmp += nRods;
            closest = int(tmp) < nRods - 1 ? int(tmp) : nRods - 1;
          } else {  // GeneralBinFinderInPhi on the PhiBorderFinder borders
            auto positive = [](double x) { return x > 0. ? x : x + 2. * M_PI; };
            for (int32_t i = 0; i < nRods && closest < 0; ++i) {
              double cur = rods[i].border(), next = rods[(i + 1) % nRods].border(), ph = phi;
              if (cur > next) {
                cur = positive(cur);
                next = positive(next);
                ph = positive(ph);
              }
              if (ph > cur - 1e-6 && ph < next)
                closest = i;
            }
            if (closest < 0)
              continue;  // legacy throws here; never seen
          }
          double fx = 0., fsx = 0.;
          bool checkNext = false;
          double dist = 0.;
          if (visitRod(closest, fx, fsx)) {
            double halfWid = rods[closest].halfX();
            if (params.btlRodsOverlap)
              halfWid *= 0.9;
            dist = fx;
            checkNext = std::abs(dist) + 3. * fsx > halfWid;
          } else {
            checkNext = true;
            double const rodPhi = std::atan2(double(rods[closest].py()), double(rods[closest].px()));
            double dphi = phi - rodPhi;
            dphi = dphi > M_PI ? dphi - 2. * M_PI : (dphi < -M_PI ? dphi + 2. * M_PI : dphi);
            dist = dphi > 0. ? -1. : 1.;
          }
          if (checkNext) {
            int32_t next = dist < 0. ? closest + 1 : closest - 1;
            next = ((next % nRods) + nRods) % nRods;
            visitRod(next, fx, fsx);
          }
        } else {
          // ================= ETL: hits in a window of cells around the crossing of the disc =================
          V3 const refPos{0., 0., params.layerRef[layer]}, refN{0., 0., 1.}, refAx{1., 0., 0.}, refAy{0., 1., 0.};
          double lx, ly, cxx, cxy, cyy, path;
          if (!localPositionAndError(st, cov, refPos, refN, refAx, refAy, lx, ly, cxx, cxy, cyy, path))
            continue;
          V3 cpos, cdir;
          moveAlong(st, path, cpos, cdir);
          // window: nSigma x (track error (+) hit error), plus the spread of the module planes around
          // the reference plane, plus 1 cm
          float const herr = params.maxHitErr[layer];
          double const dn = std::abs(cdir.z);
          double const tanIncl = std::sqrt(std::fmax(0., 1. - dn * dn)) / std::fmax(dn, 1e-3);
          double const geo = params.layerSpread[layer] * tanIncl + 1.;
          double const hx = nSig * std::sqrt(cxx + herr * herr) + geo;
          double const hy = nSig * std::sqrt(cyy + herr * herr) + geo;
          // All ETL module planes are (nearly) parallel to the disc: the Jacobian of the global (x, y)
          // crossing is computed on two z planes around the disc and interpolated in z for each module,
          // then rotated to the module axes. One crossing per module instead of eleven.
          double const sig2Ref = std::fmax(cxx, cyy) + std::abs(cxy);
          double ejx[2][5], ejy[2][5];
          double ez[2] = {0., 0.};
          bool etlFast = true;
          double const dz = std::fmax(params.layerSpread[layer], 0.5);
          for (int k = 0; k < 2; ++k) {
            ez[k] = params.layerRef[layer] + (k == 0 ? -dz : dz);
            if (!jacobianOnPlane(st, cov, V3{0., 0., ez[k]}, refN, refAx, refAy, ejx[k], ejy[k]))
              etlFast = false;
          }
          // the track on a module plane is the same for all hits of that module: keep the last one
          uint32_t lastDet = 0;
          bool lastOk = false;
          double tx = 0, ty = 0, txx = 0, txy = 0, tyy = 0, s = 0;
          auto visit = [&](int cell) {
            for (int32_t h = cells[cell].start(); h < cells[cell + 1].start(); ++h) {
              auto const hit = hits[h];
              if (hit.detId() != lastDet) {
                lastDet = hit.detId();
                lastOk = false;
                V3 const P{hit.planeX(), hit.planeY(), hit.planeZ()};
                V3 const N{hit.nX(), hit.nY(), hit.nZ()};
                V3 const AX{hit.axX(), hit.axY(), hit.axZ()};
                V3 const AY{hit.ayX(), hit.ayY(), hit.ayZ()};
                double s0;
                if (!crossPlane(st, P, N, s0))
                  continue;
                V3 p0, d0;
                moveAlong(st, s0, p0, d0);
                double const x0 = dot(sub(p0, P), AX), y0 = dot(sub(p0, P), AY);
                double const loose = std::sqrt(4. * sig2Ref);
                if (std::abs(x0) > hit.halfX() + nSig * loose || std::abs(y0) > hit.halfY() + nSig * loose)
                  continue;
                if (etlFast && std::abs(N.z) > 0.999) {
                  s = s0;
                  tx = x0;
                  ty = y0;
                  double const w = (P.z - ez[0]) / (ez[1] - ez[0]);
                  double jx[5], jy[5];
                  for (int i = 0; i < 5; ++i) {
                    double const gx = ejx[0][i] + w * (ejx[1][i] - ejx[0][i]);
                    double const gy = ejy[0][i] + w * (ejy[1][i] - ejy[0][i]);
                    jx[i] = AX.x * gx + AX.y * gy;
                    jy[i] = AY.x * gx + AY.y * gy;
                  }
                  errorFromJacobian(cov, jx, jy, txx, txy, tyy);
                } else if (!localPositionAndError(st, cov, P, N, AX, AY, tx, ty, txx, txy, tyy, s))
                  continue;
                // legacy: only modules whose bounds, widened by nSig track sigmas, contain the crossing
                if (std::abs(tx) > hit.halfX() + nSig * std::sqrt(txx) || std::abs(ty) > hit.halfY() + nSig * std::sqrt(tyy))
                  continue;
                lastOk = true;
              }
              if (lastOk)
                evalHit(h, tx, ty, txx, txy, tyy, s);
            }
          };
          int const ix0 = CellGrid::xyBin(float(cpos.x - hx)), ix1 = CellGrid::xyBin(float(cpos.x + hx));
          int const iy0 = CellGrid::xyBin(float(cpos.y - hy)), iy1 = CellGrid::xyBin(float(cpos.y + hy));
          for (int iy = iy0; iy <= iy1; ++iy)
            for (int ix = ix0; ix <= ix1; ++ix)
              visit(CellGrid::etlCell(layer, ix, iy));
        }

        // ---- 4) legacy: only the best candidate is tested against the cuts
        float const chi2Cut = layer == 0 ? params.btlChi2Cut : params.etlChi2Cut;
        float const tchi2Cut = layer == 0 ? params.btlTimeChi2Cut : params.etlTimeChi2Cut;
        if (best >= 0 && bestEst < chi2Cut && bestTchi2 < tchi2Cut) {
          records[r].hit() = best;
          records[r].est() = bestEst;
          records[r].timeChi2() = bestTchi2;
          records[r].path() = bestPath;
        }
      }
    }
  };

  // ------------------------------------------------------------------ kernel 2
  class FinishTimingKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  MtdTrackInSoA::ConstView tracks,
                                  MtdHitInSoA::ConstView hits,
                                  MtdRecordSoA::ConstView records,
                                  MtdTimingSoA::View out,
                                  MatchParams params) const {
      for (int32_t i : uniform_elements(acc, tracks.metadata().size())) {
        auto o = out[i];
        o.btlChi2() = -1.f;
        o.btlTimeChi2() = -1.f;
        o.etlChi2() = -1.f;
        o.etlTimeChi2() = -1.f;
        o.btlHit() = -1;
        o.etlHit1() = -1;
        o.etlHit2() = -1;
        o.timed() = 0;
        o.pathLength() = -1.f;
        o.tmtd() = 0.f;
        o.sigmatmtd() = -1.f;
        o.mtdPosX() = 0.f;
        o.mtdPosY() = 0.f;
        o.mtdPosZ() = 0.f;
        o.tofPi() = 0.f;
        o.tofK() = 0.f;
        o.tofP() = 0.f;
        o.sigmaTofPi() = -1.f;
        o.sigmaTofK() = -1.f;
        o.sigmaTofP() = -1.f;
        o.t0() = 0.f;
        o.sigmat0() = -1.f;
        o.beta() = 0.f;
        o.sigmaBeta() = 0.f;

        auto const t = tracks[i];
        if (!t.valid())
          continue;

        // ---- collect the matched layers
        auto const rb = records[i * kNLayers + 0];
        bool const hasBtl = rb.hit() >= 0;
        if (hasBtl) {
          o.btlHit() = hits[rb.hit()].orig();
          o.btlChi2() = rb.est();
          o.btlTimeChi2() = rb.timeChi2();
        }
        int const l1 = t.outZ() > 0.f ? 3 : 1;  // the two discs on this side
        int etl[2];
        int nEtl = 0;
        for (int l = l1; l < l1 + 2; ++l)
          if (records[i * kNLayers + l].hit() >= 0)
            etl[nEtl++] = i * kNLayers + l;
        if (nEtl > 0) {
          int bestRec = etl[0];
          if (nEtl == 2 && betterHit(records[etl[1]].est(),
                                     records[etl[1]].timeChi2(),
                                     records[etl[0]].est(),
                                     records[etl[0]].timeChi2()))
            bestRec = etl[1];
          o.etlChi2() = records[bestRec].est();
          o.etlTimeChi2() = records[bestRec].timeChi2();
          // inner = smaller |z| of the hit
          if (nEtl == 2 &&
              std::abs(hits[records[etl[0]].hit()].gz()) > std::abs(hits[records[etl[1]].hit()].gz())) {
            int const tmp = etl[0];
            etl[0] = etl[1];
            etl[1] = tmp;
          }
          o.etlHit1() = hits[records[etl[0]].hit()].orig();
          if (nEtl == 2)
            o.etlHit2() = hits[records[etl[1]].hit()].orig();
        }

        // ---- the time measurement (legacy buildTrack)
        int const nHits = (hasBtl ? 1 : 0) + nEtl;
        float thit = 0.f, terr = -1.f, gap = 0.f;
        float px = 0.f, py = 0.f, pz = 0.f;
        bool validTime = false;
        if (nHits == 1) {
          int const rec = hasBtl ? i * kNLayers : etl[0];
          auto const h = hits[records[rec].hit()];
          thit = h.time();
          terr = h.timeError();
          gap = std::abs(records[rec].path());
          px = h.gx();
          py = h.gy();
          pz = h.gz();
          validTime = true;
        } else if (nEtl == 2 && !hasBtl) {
          // combine the two ETL hits at the inner one: move the outer hit's time back along the
          // track (pion, outer momentum), compare, then average or keep the more precise one
          auto const hIn = hits[records[etl[0]].hit()];
          auto const hOut = hits[records[etl[1]].hit()];
          float const sIn = std::abs(records[etl[0]].path());
          float const sOut = std::abs(records[etl[1]].path());
          float const etlPath = std::abs(sOut - sIn);
          float const p2Out = t.outPx() * t.outPx() + t.outPy() * t.outPy() + t.outPz() * t.outPz();
          gap = sIn;
          if (etlPath > 0.f) {
            float const dtPi = etlPath * kCinv / betaOf(p2Out, kMassPi);
            float const dtP = etlPath * kCinv / betaOf(p2Out, kMassP);
            float const tBack = hOut.time() - dtPi;
            float err1 = hOut.timeError() * hOut.timeError() + (dtP - dtPi) * (dtP - dtPi);
            float err2 = hIn.timeError() * hIn.timeError();
            if (err1 > 0.f && err2 > 0.f) {
              if ((tBack - hIn.time()) * (tBack - hIn.time()) < (err1 + err2) * params.etlTimeChi2Cut) {
                float const w1 = 1.f / err1, w2 = 1.f / err2;
                float const e2 = 1.f / (w1 + w2);
                thit = (tBack * w1 + hIn.time() * w2) * e2;
                terr = std::sqrt(e2);
                px = hIn.gx();
                py = hIn.gy();
                pz = hIn.gz();
              } else if (err1 <= err2) {  // legacy: position left at (0,0,0) in these two branches
                thit = tBack;
                terr = std::sqrt(err1);
              } else {
                thit = hIn.time();
                terr = hIn.timeError();
              }
              validTime = true;
            }
          }
        }
        if (!validTime)
          continue;  // BTL+ETL or 3 hits: legacy "anomalous pattern", no time

        TofInput const tin = tofInputOf(t, gap);
        float const tofPi = tofFor(tin, kMassPi), tofK = tofFor(tin, kMassK), tofP = tofFor(tin, kMassP);
        float const p2 = tin.p2Pca;
        o.timed() = 1;
        o.pathLength() = tin.trackerPath + tin.gapPath;
        o.tmtd() = thit;
        o.sigmatmtd() = terr;
        o.mtdPosX() = px;
        o.mtdPosY() = py;
        o.mtdPosZ() = pz;
        o.tofPi() = tofPi;
        o.tofK() = tofK;
        o.tofP() = tofP;
        o.sigmaTofPi() = sigmaTofFor(tin, kMassPi);
        o.sigmaTofK() = sigmaTofFor(tin, kMassK);
        o.sigmaTofP() = sigmaTofFor(tin, kMassP);
        o.t0() = thit - tofPi;
        o.sigmat0() = std::sqrt(terr * terr + (tofP - tofPi) * (tofP - tofPi));
        o.beta() = betaOf(p2, kMassPi);
        o.sigmaBeta() = betaOf(p2, kMassP) - betaOf(p2, kMassPi);
      }
    }
  };

  // ------------------------------------------------------------------ kernel 3
  // MTDTrackQualityMVA: the same cuts, inputs and default (-1)
  class QualityKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  MtdTrackInSoA::ConstView tracks,
                                  MtdBdtSoA::ConstView bdt,
                                  MtdTimingSoA::View out) const {
      for (int32_t i : uniform_elements(acc, tracks.metadata().size())) {
        auto const t = tracks[i];
        auto o = out[i];
        o.mva() = -1.f;
        double const aeta = std::abs(t.eta());
        if (aeta > 3. || (aeta <= 1.5 && t.pt() < 0.7) || (aeta > 1.5 && t.pt() < 0.2))
          continue;
        if (!(o.tmtd() > 0))
          continue;
        // outermost hit position, as filled by MtdTimingToValueMaps
        bool const btl = o.btlHit() >= 0, etl = o.etlHit1() >= 0;
        float outermost;
        if (btl || etl)
          outermost = btl ? t.outerRadius() : t.outerZ();
        else
          outermost = aeta < 1.5f ? t.outerRadius() : t.outerZ();
        float const vars[kNBdtVars] = {float(t.pt()),
                                       float(t.eta()),
                                       t.phi(),
                                       t.dzBS(),
                                       t.dxyBS(),
                                       t.chi2(),
                                       t.ndof(),
                                       float(t.npixBarrel()),
                                       float(t.npixEndcap()),
                                       o.btlChi2(),
                                       o.btlTimeChi2(),
                                       o.etlChi2(),
                                       o.etlTimeChi2(),
                                       o.tmtd(),
                                       o.sigmatmtd(),
                                       o.pathLength(),
                                       outermost};
        o.mva() = bdtQuality(bdt, vars);
      }
    }
  };

  void runMtdTiming(Queue& queue,
                    MatchParams const& params,
                    MtdEventInputSoA::ConstView input,
                    MtdBtlGeometrySoA::ConstView geometry,
                    MtdBdtSoA::ConstView bdt,
                    RecordDevice& records,
                    MtdTimingDevice& out) {
    auto const tracks = input.tracks();
    auto const hits = input.hits();
    int32_t const nTracks = tracks.metadata().size();
    if (nTracks == 0)
      return;
    constexpr uint32_t threads = 64;
    auto const wd1 = make_workdiv<Acc1D>(divide_up_by(nTracks * kNLayers, threads), threads);
    alpaka::exec<Acc1D>(
        queue,
        wd1,
        MatchHitsKernel{},
        tracks,
        hits,
        input.cells(),
        geometry.rods(),
        geometry.modules(),
        records.view(),
        params);
    auto const wd2 = make_workdiv<Acc1D>(divide_up_by(nTracks, threads), threads);
    alpaka::exec<Acc1D>(queue,
                        wd2,
                        FinishTimingKernel{},
                        tracks,
                        hits,
                        records.const_view(),
                        out.view(),
                        params);
    alpaka::exec<Acc1D>(queue, wd2, QualityKernel{}, tracks, bdt, out.view());
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming
