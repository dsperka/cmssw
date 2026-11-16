#ifndef RecoMTD_TrackExtenderAlpaka_interface_MtdTimingInputSoA_h
#define RecoMTD_TrackExtenderAlpaka_interface_MtdTimingInputSoA_h

// Internal tables of the MTD timing module (not event products).
//   Tracks: one row per input track.
//   Hits:   one row per MTD hit, sorted by (layer, cell); each hit carries its module plane.
//   Cells:  start index of every cell in the hit table (CSR), size nCells + 1.
//   Records: scratch, one row per (track, layer) with the best hit found in that layer.

#include <cmath>
#include <cstdint>

#include <Eigen/Core>

#include <alpaka/core/Common.hpp>

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace mtdtiming {

  using Vector15f = Eigen::Matrix<float, 15, 1>;

  GENERATE_SOA_LAYOUT(MtdTrackInLayout,
                      SOA_COLUMN(int8_t, valid),  // 0 if the track has no TrackExtra (no outer state)
                      SOA_COLUMN(float, charge),
                      SOA_COLUMN(float, bz),  // field at the outer state (T)
                      SOA_COLUMN(float, outX),
                      SOA_COLUMN(float, outY),
                      SOA_COLUMN(float, outZ),
                      SOA_COLUMN(float, outPx),
                      SOA_COLUMN(float, outPy),
                      SOA_COLUMN(float, outPz),
                      SOA_EIGEN_COLUMN(Vector15f, outCov),  // curvilinear 5x5, lower triangle
                      SOA_COLUMN(float, x0),                // reference point (PCA)
                      SOA_COLUMN(float, y0),
                      SOA_COLUMN(float, z0),
                      SOA_COLUMN(float, ptPca),
                      SOA_COLUMN(float, pPca),
                      SOA_COLUMN(float, cosThetaPca),
                      SOA_COLUMN(float, sigmaPPca))

  GENERATE_SOA_LAYOUT(MtdHitInLayout,
                      SOA_COLUMN(int8_t, layer),  // 0 = BTL, 1/2 = ETL -z disc 1/2, 3/4 = ETL +z disc 1/2
                      SOA_COLUMN(uint32_t, detId),
                      SOA_COLUMN(int32_t, orig),  // index of the hit in the input MTDTrackingDetSetVector data
                      SOA_COLUMN(float, planeX),  // module plane: position and local axes
                      SOA_COLUMN(float, planeY),
                      SOA_COLUMN(float, planeZ),
                      SOA_COLUMN(float, axX),
                      SOA_COLUMN(float, axY),
                      SOA_COLUMN(float, axZ),
                      SOA_COLUMN(float, ayX),
                      SOA_COLUMN(float, ayY),
                      SOA_COLUMN(float, ayZ),
                      SOA_COLUMN(float, nX),
                      SOA_COLUMN(float, nY),
                      SOA_COLUMN(float, nZ),
                      SOA_COLUMN(float, halfX),  // module half-width / half-length (legacy compatibility test)
                      SOA_COLUMN(float, halfY),
                      SOA_COLUMN(float, lx),  // hit local position and error
                      SOA_COLUMN(float, ly),
                      SOA_COLUMN(float, exx),
                      SOA_COLUMN(float, exy),
                      SOA_COLUMN(float, eyy),
                      SOA_COLUMN(float, time),
                      SOA_COLUMN(float, timeError),
                      SOA_COLUMN(float, gx),  // hit global position
                      SOA_COLUMN(float, gy),
                      SOA_COLUMN(float, gz))

  GENERATE_SOA_LAYOUT(MtdCellLayout, SOA_COLUMN(int32_t, start))

  // BTL geometry as the legacy navigation sees it (MTDTrayBarrelLayer): "rods" (rows of modules
  // along z, sorted in phi) and their modules (in rod order). For BTL the hit "cell" is the module
  // index, so the cell index gives the hits of each module.
  GENERATE_SOA_LAYOUT(MtdPlaneLayout,
                      SOA_COLUMN(float, px),
                      SOA_COLUMN(float, py),
                      SOA_COLUMN(float, pz),
                      SOA_COLUMN(float, axX),
                      SOA_COLUMN(float, axY),
                      SOA_COLUMN(float, axZ),
                      SOA_COLUMN(float, ayX),
                      SOA_COLUMN(float, ayY),
                      SOA_COLUMN(float, ayZ),
                      SOA_COLUMN(float, nX),
                      SOA_COLUMN(float, nY),
                      SOA_COLUMN(float, nZ),
                      SOA_COLUMN(float, halfX),
                      SOA_COLUMN(float, halfY),
                      SOA_COLUMN(int32_t, first),  // rods: first module; modules: rod index
                      SOA_COLUMN(int32_t, count),  // rods: number of modules; modules: unused
                      SOA_COLUMN(float, border))   // rods: lower phi border (PhiBorderFinder); modules: unused

  GENERATE_SOA_LAYOUT(MtdRecordLayout,
                      SOA_COLUMN(int32_t, hit),  // best hit passing the cuts, -1 if none
                      SOA_COLUMN(float, est),
                      SOA_COLUMN(float, timeChi2),
                      SOA_COLUMN(float, path))  // signed path from the outer state to the hit plane

  using MtdTrackInSoA = MtdTrackInLayout<>;
  using MtdHitInSoA = MtdHitInLayout<>;
  using MtdCellSoA = MtdCellLayout<>;
  using MtdRecordSoA = MtdRecordLayout<>;
  using MtdPlaneSoA = MtdPlaneLayout<>;

  constexpr int kNLayers = 5;

  // ---- the cell grid used to find nearby hits (about 2 cm cells)
  // BTL: (phi, z) cells. ETL: (x, y) cells, one grid per disc.
  struct CellGrid {
    static constexpr int kBtlNPhi = 360;
    static constexpr int kBtlNZ = 300;
    static constexpr float kBtlZMin = -300.f, kBtlZMax = 300.f;
    static constexpr int kEtlN = 150;
    static constexpr float kEtlMin = -150.f, kEtlMax = 150.f;
    static constexpr int kBtlCells = kBtlNPhi * kBtlNZ;
    static constexpr int kEtlCells = kEtlN * kEtlN;
    static constexpr int kNCells = kBtlCells + 4 * kEtlCells;

    ALPAKA_FN_HOST_ACC static int clampBin(float v, float lo, float hi, int n) {
      int b = int((v - lo) / (hi - lo) * n);
      return b < 0 ? 0 : (b >= n ? n - 1 : b);
    }
    ALPAKA_FN_HOST_ACC static int phiBin(float phi) {  // phi in [-pi, pi]
      return clampBin(phi, -float(M_PI), float(M_PI), kBtlNPhi);
    }
    ALPAKA_FN_HOST_ACC static int zBin(float z) { return clampBin(z, kBtlZMin, kBtlZMax, kBtlNZ); }
    ALPAKA_FN_HOST_ACC static int xyBin(float v) { return clampBin(v, kEtlMin, kEtlMax, kEtlN); }
    ALPAKA_FN_HOST_ACC static int btlCell(int iphi, int iz) { return iz * kBtlNPhi + iphi; }
    ALPAKA_FN_HOST_ACC static int etlCell(int layer, int ix, int iy) {
      return kBtlCells + (layer - 1) * kEtlCells + iy * kEtlN + ix;
    }
    ALPAKA_FN_HOST_ACC static int cellOf(int layer, float x, float y, float z) {
      if (layer == 0)
        return btlCell(phiBin(std::atan2(y, x)), zBin(z));
      return etlCell(layer, xyBin(x), xyBin(y));
    }
  };

  // ---- per-event parameters passed to the kernels by value
  struct MatchParams {
    float btlChi2Cut, btlTimeChi2Cut, etlChi2Cut, etlTimeChi2Cut;
    float estimatorMaxChi2;  // 500: candidates above this are ignored (legacy estimator)
    float nSigmaWindow;      // half-size of the search window in track sigmas (legacy: 10)
    float bsTimeSpread;      // 0.2 ns
    int debug;
    int skipCandidates;      // validation/timing only               // >0: print the per-layer candidate summary (validation only)
    // per layer, measured from the hits of the event on the host
    float layerRef[kNLayers];     // BTL: mean module radius; ETL: mean module z
    float layerSpread[kNLayers];  // max |r - ref| (BTL) or |z - ref| (ETL) of the modules
    float maxHitErr[kNLayers];    // max sqrt(hit local error) in the layer
    float btlHalfLength;          // BTL layer cylinder: half-length and half-thickness
    float btlHalfThickness;
    int btlNRods;                 // rod bin finder (PeriodicBinFinderInPhi)
    float btlRodPhiOffset, btlRodPhiStep;
    int btlRodsOverlap;           // PhiBorderFinder::isPhiOverlapping
    int btlRodsPeriodic;          // PhiBorderFinder::isPhiPeriodic (else the general bin finder on the borders)
  };

}  // namespace mtdtiming

#endif
