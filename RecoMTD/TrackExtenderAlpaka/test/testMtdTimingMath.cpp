// Step 2 unit test: compare the simple helix math (MtdTimingMath.h) with the legacy
// AnalyticalPropagator, which is the geometric part of PropagatorWithMaterialForMTD used by the
// legacy hit matching (no material effects exist on MTD surfaces).
//
// For random track states near the tracker outer surface and random MTD-like planes
// (BTL-like tilted planes at R ~ 116 cm, ETL-like planes at |z| ~ 300 cm) it compares:
//   path length, local x/y of the crossing, and the 2x2 local position error.
// Exit code 0 if all agree within tolerances.

#include <cmath>
#include <cstdio>
#include <random>

#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToLocal.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "RecoMTD/TrackExtenderAlpaka/interface/MtdTimingMath.h"

namespace {
  class UniformField : public MagneticField {
  public:
    explicit UniformField(float bz) : bz_(bz) { setNominalValue(); }
    GlobalVector inTesla(const GlobalPoint&) const override { return GlobalVector(0, 0, bz_); }
    bool isDefined(const GlobalPoint&) const override { return true; }
    int nominalValue() const { return 38; }

  private:
    float bz_;
  };
}  // namespace

int main() {
  using namespace mtdtiming;
  constexpr float bz = 3.8112f;
  UniformField field(bz);
  AnalyticalPropagator prop(&field, anyDirection, 1.6);

  std::mt19937 rng(12345);
  std::uniform_real_distribution<double> uni(0., 1.);
  auto U = [&](double a, double b) { return a + (b - a) * uni(rng); };

  int nTried = 0, nBoth = 0, nOnlyOne = 0, nBad = 0;
  double maxDPath = 0, maxDPos = 0, maxRelErr = 0;
  for (int i = 0; i < 20000; ++i) {
    bool const barrel = i % 2 == 0;
    // ---- random state at the tracker outer surface
    double const pt = std::exp(U(std::log(0.6), std::log(20.)));
    double const q = uni(rng) < 0.5 ? -1. : 1.;
    double phi = U(-M_PI, M_PI);
    double eta, x, y, z;
    if (barrel) {
      eta = U(-1.4, 1.4);
      double const r = U(100., 110.);
      x = r * std::cos(phi + U(-0.2, 0.2));
      y = r * std::sin(phi + U(-0.2, 0.2));
      z = r * std::sinh(eta);
    } else {
      eta = (uni(rng) < 0.5 ? -1 : 1) * U(1.7, 2.9);
      double const zz = (eta > 0 ? 1 : -1) * U(270., 285.);
      double const r = std::abs(zz / std::sinh(eta));
      x = r * std::cos(phi);
      y = r * std::sin(phi);
      z = zz;
    }
    double const pz = pt * std::sinh(eta);
    GlobalPoint gp(x, y, z);
    GlobalVector gm(pt * std::cos(phi), pt * std::sin(phi), pz);
    GlobalTrajectoryParameters gtp(gp, gm, int(q), &field);
    // random positive-definite curvilinear covariance (typical outer-state errors)
    AlgebraicSymMatrix55 c;
    double const s[5] = {0.02 / gm.mag(), 2e-3, 2e-3, 0.05, 0.08};
    for (int a = 0; a < 5; ++a)
      c(a, a) = s[a] * s[a];
    c(1, 4) = 0.3 * s[1] * s[4];
    c(2, 3) = -0.4 * s[2] * s[3];
    c(0, 2) = 0.2 * s[0] * s[2];
    FreeTrajectoryState fts(gtp, CurvilinearTrajectoryError(c));

    // ---- random plane
    GlobalPoint ppos;
    GlobalVector ax, ay;
    if (barrel) {
      double const r = U(114.5, 117.);
      double const ph = std::atan2(y, x) + U(-0.05, 0.05);
      double const tilt = U(-0.3, 0.3);  // rotation of the module around z
      GlobalVector radial(std::cos(ph + tilt), std::sin(ph + tilt), 0.);
      ppos = GlobalPoint(r * std::cos(ph), r * std::sin(ph), z * r / std::hypot(x, y) + U(-5, 5));
      ax = GlobalVector(-radial.y(), radial.x(), 0.);
      ay = GlobalVector(0, 0, 1);
    } else {
      double const zp = (z > 0 ? 1 : -1) * U(298., 304.);
      ppos = GlobalPoint(x * zp / z + U(-3, 3), y * zp / z + U(-3, 3), zp);
      double const a = U(-M_PI, M_PI);
      ax = GlobalVector(std::cos(a), std::sin(a), 0.);
      ay = GlobalVector(-std::sin(a), std::cos(a), 0.);
    }
    Surface::RotationType rot(ax.basicVector(), ay.basicVector());
    auto plane = Plane::build(ppos, rot);

    // ---- legacy
    ++nTried;
    Propagator const& baseProp = prop;
    auto res = baseProp.propagateWithPath(fts, *plane);
    bool const okLegacy = res.first.isValid();

    // ---- new
    State st{{x, y, z}, {gm.x() / gm.mag(), gm.y() / gm.mag(), gm.z() / gm.mag()}, gm.mag(), q, bz};
    float cov15[15];
    int k = 0;
    for (int a = 0; a < 5; ++a)
      for (int b = 0; b <= a; ++b)
        cov15[k++] = c(a, b);
    V3 const P{ppos.x(), ppos.y(), ppos.z()};
    V3 const AX{rot.xx(), rot.xy(), rot.xz()}, AY{rot.yx(), rot.yy(), rot.yz()}, AZ{rot.zx(), rot.zy(), rot.zz()};
    double lx, ly, cxx, cxy, cyy, path;
    bool const okNew = localPositionAndError(st, cov15, P, AZ, AX, AY, lx, ly, cxx, cxy, cyy, path);

    if (okLegacy != okNew) {
      ++nOnlyOne;
      continue;
    }
    if (!okLegacy)
      continue;
    ++nBoth;
    auto const& t = res.first;
    double const dPath = std::abs(path - res.second);
    double const dPos = std::hypot(lx - t.localPosition().x(), ly - t.localPosition().y());
    auto const e = t.localError().positionError();
    double const rel = std::max({std::abs(cxx - e.xx()) / e.xx(),
                                 std::abs(cyy - e.yy()) / e.yy(),
                                 std::abs(cxy - e.xy()) / std::sqrt(e.xx() * e.yy())});
    maxDPath = std::max(maxDPath, dPath);
    maxDPos = std::max(maxDPos, dPos);
    maxRelErr = std::max(maxRelErr, rel);
    if (dPath > 1e-3 || dPos > 1e-3 || rel > 1e-3) {
      if (nBad < 10)
        std::printf("MISMATCH %s pt=%.2f q=%+.0f path %.5f vs %.5f  pos (%.5f,%.5f) vs (%.5f,%.5f)  "
                    "err xx %.4g vs %.4g yy %.4g vs %.4g xy %.4g vs %.4g\n",
                    barrel ? "BTL" : "ETL", pt, q, path, res.second, lx, ly, t.localPosition().x(),
                    t.localPosition().y(), cxx, e.xx(), cyy, e.yy(), cxy, e.xy());
      ++nBad;
    }
  }
  std::printf("tried %d, both found %d, only one found %d, mismatches %d\n", nTried, nBoth, nOnlyOne, nBad);
  std::printf("max |d path| %.2e cm, max |d pos| %.2e cm, max rel. error difference %.2e\n", maxDPath, maxDPos,
              maxRelErr);
  bool const pass = nBad == 0 && nOnlyOne <= nTried / 1000;
  std::printf("%s\n", pass ? "PASS" : "FAIL");
  return pass ? 0 : 1;
}
