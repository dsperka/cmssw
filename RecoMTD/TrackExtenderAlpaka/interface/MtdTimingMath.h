#ifndef RecoMTD_TrackExtenderAlpaka_interface_MtdTimingMath_h
#define RecoMTD_TrackExtenderAlpaka_interface_MtdTimingMath_h

// Physics math for the MTD track timing, shared by the host and the device code.
// Plain functions, no CMSSW or Alpaka types except the ALPAKA_FN_HOST_ACC qualifier.
//
// Units: cm, GeV, Tesla, ns.
//
// What is here:
//   * Helix: a charged track in a uniform field Bz (the field at the starting point, like the
//     legacy AnalyticalPropagator). Position and direction after a path length s.
//   * Crossing of the helix with a plane (any orientation; covers tilted modules).
//   * The 2x2 error of the track position on a plane: we move each of the 5 track parameters a
//     little, propagate again, and take differences (numerical Jacobian). Simple and easy to check.
//   * The tracker path length estimate without a Trajectory (rule "C" from the step 0 study).
//   * Time of flight for pi/K/p and the legacy time chi2 rule (including the "heavy" choice).

#include <cmath>

#include <alpaka/core/Common.hpp>

namespace mtdtiming {

  constexpr double kBFieldToCurvature = 0.00299792458;  // GeV / (T cm)
  constexpr float kCcmNs = 29.9792458f;                 // speed of light, cm/ns
  constexpr float kCinv = 1.f / kCcmNs;
  constexpr float kMassPi = 0.13957018f;
  constexpr float kMassK = 0.493677f;
  constexpr float kMassP = 0.9382720813f;

  // ---------------------------------------------------------------- small vector helpers
  struct V3 {
    double x, y, z;
  };
  ALPAKA_FN_HOST_ACC inline V3 add(V3 a, V3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
  ALPAKA_FN_HOST_ACC inline V3 sub(V3 a, V3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
  ALPAKA_FN_HOST_ACC inline V3 mul(V3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
  ALPAKA_FN_HOST_ACC inline double dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

  // ---------------------------------------------------------------- the track state
  // Position, unit direction, |p|, charge and the field Bz used for the helix.
  struct State {
    V3 pos;
    V3 dir;  // unit vector along the momentum
    double p;
    double q;
    double bz;
  };

  // Turning rate of the transverse direction per unit path length (rad/cm).
  ALPAKA_FN_HOST_ACC inline double omega(State const& s) { return -s.q * kBFieldToCurvature * s.bz / s.p; }

  // Position and direction after path length s along the helix (s may be negative).
  ALPAKA_FN_HOST_ACC inline void moveAlong(State const& st, double s, V3& pos, V3& dir) {
    double const w = omega(st);
    double const a = w * s;
    double sa, ca;
    double const tx = st.dir.x, ty = st.dir.y, tz = st.dir.z;
    // sin(a)/w and (1-cos(a))/w, with a series for small a to avoid 0/0
    double sinOverW, oneMinusCosOverW;
    if (std::abs(a) < 1e-4) {
      sinOverW = s * (1. - a * a / 6.);
      oneMinusCosOverW = s * a * (0.5 - a * a / 24.);
    } else {
      sinOverW = std::sin(a) / w;
      oneMinusCosOverW = (1. - std::cos(a)) / w;
    }
    sa = std::sin(a);
    ca = std::cos(a);
    pos.x = st.pos.x + tx * sinOverW - ty * oneMinusCosOverW;
    pos.y = st.pos.y + ty * sinOverW + tx * oneMinusCosOverW;
    pos.z = st.pos.z + tz * s;
    dir.x = tx * ca - ty * sa;
    dir.y = ty * ca + tx * sa;
    dir.z = tz;
  }

  // ---------------------------------------------------------------- crossing with a plane
  // Plane: a point on it and its unit normal. Returns false if no crossing was found.
  // Starts from the straight-line solution and refines with Newton steps (the nearest crossing,
  // in either direction, like the legacy "anyDirection" propagator). Rejects turns above maxDPhi.
  ALPAKA_FN_HOST_ACC inline bool crossPlane(State const& st, V3 planePos, V3 normal, double& s, double maxDPhi = 1.6) {
    double const d0 = dot(normal, st.dir);
    if (std::abs(d0) < 1e-9)
      return false;
    s = -dot(normal, sub(st.pos, planePos)) / d0;
    V3 pos, dir;
    for (int it = 0; it < 20; ++it) {
      moveAlong(st, s, pos, dir);
      double const f = dot(normal, sub(pos, planePos));
      double const df = dot(normal, dir);
      if (std::abs(df) < 1e-9)
        return false;
      double const ds = f / df;
      s -= ds;
      if (std::abs(ds) < 1e-6)
        return std::abs(omega(st) * s) <= maxDPhi;
    }
    return false;
  }

  // Crossing with a cylinder of radius r around the z axis (BTL), forward or backward,
  // nearest to the start. Used only to centre the hit search window.
  ALPAKA_FN_HOST_ACC inline bool crossCylinder(State const& st, double r, double& s, double maxDPhi = 1.6) {
    // straight-line start: |pos + s dir|_T = r
    double const a = st.dir.x * st.dir.x + st.dir.y * st.dir.y;
    if (a < 1e-12)
      return false;
    double const b = 2. * (st.pos.x * st.dir.x + st.pos.y * st.dir.y);
    double const c = st.pos.x * st.pos.x + st.pos.y * st.pos.y - r * r;
    double const disc = b * b - 4. * a * c;
    if (disc < 0.)
      return false;
    double const s1 = (-b + std::sqrt(disc)) / (2. * a);
    double const s2 = (-b - std::sqrt(disc)) / (2. * a);
    s = std::abs(s1) < std::abs(s2) ? s1 : s2;
    V3 pos, dir;
    for (int it = 0; it < 20; ++it) {
      moveAlong(st, s, pos, dir);
      double const f = pos.x * pos.x + pos.y * pos.y - r * r;
      double const df = 2. * (pos.x * dir.x + pos.y * dir.y);
      if (std::abs(df) < 1e-9)
        return false;
      double const ds = f / df;
      s -= ds;
      if (std::abs(ds) < 1e-6)
        return std::abs(omega(st) * s) <= maxDPhi;
    }
    return false;
  }

  // ---------------------------------------------------------------- track errors on a plane
  // Curvilinear parameters (CMSSW convention): (q/p, lambda, phi, xT, yT), with the position offset
  // xT * U + yT * V, U = (-sin phi, cos phi, 0), V = T x U.
  // cov15 = the 5x5 symmetric curvilinear covariance, lower triangle row by row:
  // (00, 10, 11, 20, 21, 22, 30, 31, 32, 33, 40, 41, 42, 43, 44).
  ALPAKA_FN_HOST_ACC inline double covAt(float const* cov15, int i, int j) {
    if (i < j) {
      int t = i;
      i = j;
      j = t;
    }
    return cov15[i * (i + 1) / 2 + j];
  }

  // Build a state from curvilinear parameters around a reference state.
  ALPAKA_FN_HOST_ACC inline State shifted(State const& ref, double const d[5]) {
    double const qop0 = ref.q / ref.p;
    double const lam0 = std::asin(ref.dir.z);
    double const phi0 = std::atan2(ref.dir.y, ref.dir.x);
    double const qop = qop0 + d[0];
    double const lam = lam0 + d[1];
    double const phi = phi0 + d[2];
    State s = ref;
    s.p = std::abs(ref.q / qop);
    s.dir = {std::cos(lam) * std::cos(phi), std::cos(lam) * std::sin(phi), std::sin(lam)};
    // frame of the reference state
    double const sp = std::sin(phi0), cp = std::cos(phi0), sl = std::sin(lam0), cl = std::cos(lam0);
    V3 const u{-sp, cp, 0.};
    V3 const v{-sl * cp, -sl * sp, cl};  // T x U
    s.pos = add(ref.pos, add(mul(u, d[3]), mul(v, d[4])));
    return s;
  }

  // Derivatives (2x5) of the local (x, y) crossing point on a plane with respect to the 5
  // curvilinear track parameters, by central differences (each parameter moved by 1e-3 of its error).
  ALPAKA_FN_HOST_ACC inline bool jacobianOnPlane(
      State const& st, float const* cov15, V3 planePos, V3 normal, V3 ax, V3 ay, double jx[5], double jy[5]) {
    V3 pos, dir;
    for (int i = 0; i < 5; ++i) {
      double const sig = std::sqrt(covAt(cov15, i, i));
      double const h = sig > 0. ? 1e-3 * sig : 1e-9;
      double d[5] = {0., 0., 0., 0., 0.};
      double sx;
      d[i] = h;
      State const sp = shifted(st, d);
      if (!crossPlane(sp, planePos, normal, sx))
        return false;
      moveAlong(sp, sx, pos, dir);
      double const xp = dot(sub(pos, planePos), ax), yp = dot(sub(pos, planePos), ay);
      d[i] = -h;
      State const sm = shifted(st, d);
      if (!crossPlane(sm, planePos, normal, sx))
        return false;
      moveAlong(sm, sx, pos, dir);
      double const xm = dot(sub(pos, planePos), ax), ym = dot(sub(pos, planePos), ay);
      jx[i] = (xp - xm) / (2. * h);
      jy[i] = (yp - ym) / (2. * h);
    }
    return true;
  }

  // 2x2 error J C J^T.
  ALPAKA_FN_HOST_ACC inline void errorFromJacobian(
      float const* cov15, double const jx[5], double const jy[5], double& cxx, double& cxy, double& cyy) {
    cxx = cxy = cyy = 0.;
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j < 5; ++j) {
        double const c = covAt(cov15, i, j);
        cxx += jx[i] * c * jx[j];
        cxy += jx[i] * c * jy[j];
        cyy += jy[i] * c * jy[j];
      }
  }

  // Track position on a plane in the plane's local (x, y) coordinates, and its 2x2 error.
  // ax, ay = the plane's local x and y axes (unit vectors). Returns false if the crossing fails.
  ALPAKA_FN_HOST_ACC inline bool localPositionAndError(State const& st,
                                                       float const* cov15,
                                                       V3 planePos,
                                                       V3 normal,
                                                       V3 ax,
                                                       V3 ay,
                                                       double& lx,
                                                       double& ly,
                                                       double& cxx,
                                                       double& cxy,
                                                       double& cyy,
                                                       double& path) {
    double s;
    if (!crossPlane(st, planePos, normal, s))
      return false;
    V3 pos, dir;
    moveAlong(st, s, pos, dir);
    lx = dot(sub(pos, planePos), ax);
    ly = dot(sub(pos, planePos), ay);
    path = s;
    double jx[5], jy[5];
    if (!jacobianOnPlane(st, cov15, planePos, normal, ax, ay, jx, jy))
      return false;
    errorFromJacobian(cov15, jx, jy, cxx, cxy, cyy);
    return true;
  }

  // Legacy Chi2MeasurementEstimator for a 2D hit: r^T (V_hit + C_track)^-1 r.
  ALPAKA_FN_HOST_ACC inline double chi2Local(
      double dx, double dy, double cxx, double cxy, double cyy) {
    double const det = cxx * cyy - cxy * cxy;
    if (!(det > 0.))
      return 1e30;
    return (dx * dx * cyy - 2. * dx * dy * cxy + dy * dy * cxx) / det;
  }

  // ---------------------------------------------------------------- tracker path (no Trajectory)
  // Path length from the point of closest approach to the outermost tracker state.
  // Central tracks: transverse helix arc between the two points, with the average pt and angle.
  // Forward tracks: delta z over the average |cos(theta)|. (Step 0: equal to the exact
  // Trajectory path in BTL, about 1 ps worse in ETL.)
  ALPAKA_FN_HOST_ACC inline double trackerPath(double x0,
                                               double y0,
                                               double z0,
                                               double ptPca,
                                               double cosThetaPca,  // signed pz/p at the PCA
                                               V3 outPos,
                                               V3 outMom,
                                               double bz) {
    double const pOut = std::sqrt(dot(outMom, outMom));
    double const ptOut = std::sqrt(outMom.x * outMom.x + outMom.y * outMom.y);
    double const cAvg = 0.5 * (std::abs(cosThetaPca) + std::abs(outMom.z) / pOut);
    if (cAvg >= 0.5)
      return std::abs(outPos.z - z0) / cAvg;
    double const ptAvg = 0.5 * (ptPca + ptOut);
    double const radius = ptAvg / (kBFieldToCurvature * std::abs(bz));
    double const chord = std::sqrt((outPos.x - x0) * (outPos.x - x0) + (outPos.y - y0) * (outPos.y - y0));
    double const ratio = chord / (2. * radius);
    double const arcT = 2. * radius * std::asin(ratio < 1. ? ratio : 1.);
    double const sinAvg = std::sqrt(1. - cAvg * cAvg);
    return arcT / sinAvg;
  }

  // ---------------------------------------------------------------- time of flight
  ALPAKA_FN_HOST_ACC inline float betaOf(float p2, float m) {
    float const gammasq = 1.f + p2 / (m * m);
    return std::sqrt(1.f - 1.f / gammasq);
  }

  // Inputs for the TOF of one track: the tracker part uses two momenta (PCA and outer) with
  // weight w on the PCA one (w = 0.37, calibrated on the legacy segment sum); the gap from the
  // outer tracker state to the MTD hit uses the outer momentum.
  struct TofInput {
    float trackerPath;  // cm
    float gapPath;      // cm
    float p2Pca, p2Out;
    float sigmaPPca;  // GeV, error of |p| at the PCA
  };
  constexpr float kTofWeightPca = 0.37f;

  ALPAKA_FN_HOST_ACC inline float tofFor(TofInput const& in, float m) {
    float const bPca = betaOf(in.p2Pca, m), bOut = betaOf(in.p2Out, m);
    return in.trackerPath * kCinv * (kTofWeightPca / bPca + (1.f - kTofWeightPca) / bOut) +
           in.gapPath * kCinv / bOut;
  }

  // Error of the TOF from the momentum error (fully correlated along the track, legacy kCost).
  ALPAKA_FN_HOST_ACC inline float sigmaTofFor(TofInput const& in, float m) {
    float const p2 = in.p2Pca;
    return (in.trackerPath + in.gapPath) * kCinv * in.sigmaPPca * m * m / (p2 * std::sqrt(p2 + m * m));
  }

  // Legacy time chi2 used to rank candidate hits (computeTrackTofPidInfo with addPIDError=false):
  // pion hypothesis, unless the hit is much more compatible with a kaon or proton ("heavy").
  // tVtxErr = beam spot time spread (0.2 ns). Computed in double so exp() cannot underflow to NaN.
  ALPAKA_FN_HOST_ACC inline float timeChi2(TofInput const& in, float tHit, float tHitErr, float tVtxErr) {
    double const base = double(tHitErr) * tHitErr + double(tVtxErr) * tVtxErr;
    double chi2[3];
    float const masses[3] = {kMassPi, kMassK, kMassP};
    for (int i = 0; i < 3; ++i) {
      double const dt = double(tHit) - tofFor(in, masses[i]);
      double const st = sigmaTofFor(in, masses[i]);
      chi2[i] = dt * dt / (base + st * st);
    }
    // probabilities relative to the smallest chi2 (same ratios as the legacy, no underflow)
    double const cmin = std::fmin(chi2[0], std::fmin(chi2[1], chi2[2]));
    double const rpi = std::exp(-0.5 * (chi2[0] - cmin));
    double const rk = std::exp(-0.5 * (chi2[1] - cmin));
    double const rp = std::exp(-0.5 * (chi2[2] - cmin));
    double const probHeavy = 1. - rpi / (rpi + rk + rp);
    if (probHeavy > 0.75)
      return float(chi2[1] < chi2[2] ? chi2[1] : chi2[2]);
    return float(chi2[0]);
  }

  // Legacy ordering of two candidate hits (MTDHitMatchingInfo::operator<).
  ALPAKA_FN_HOST_ACC inline bool betterHit(float est1, float tchi1, float est2, float tchi2) {
    if (tchi1 < 10.f && tchi2 < 10.f)
      return est1 + 3.f * tchi1 < est2 + 3.f * tchi2;
    return est1 + 8.f * tchi1 < est2 + 8.f * tchi2;
  }

}  // namespace mtdtiming

#endif
