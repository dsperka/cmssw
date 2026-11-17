#ifndef DataFormats_MtdTimingSoA_interface_MtdEventInputSoA_h
#define DataFormats_MtdTimingSoA_interface_MtdEventInputSoA_h

// Everything the device chain (MTD timing, then vertices) needs from one event, in one buffer, so
// that it goes to the device with a single copy. Built on the host by MtdEventInputProducer.
//
//   tracks    one row per input reco::Track (same order): the state at the outermost hit, for MTD
//   hits      one row per MTD hit, sorted by search cell; each hit carries its module plane
//   cells     start of every search cell in the hit table, size nCells + 1
//   layers    one row per MTD layer (BTL, ETL -z disc 1/2, ETL +z disc 1/2), measured on the hits
//   vtxTracks the tracks selected for vertexing, sorted by z at the beam line (TrackForVertexSoA)

#include <cstdint>

#include <Eigen/Core>

#include "DataFormats/SoATemplate/interface/SoABlocks.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/VertexSoA/interface/TrackForVertexSoA.h"

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
                      SOA_COLUMN(float, sigmaPPca),
                      // inputs of the MTD track-quality BDT (MTDTrackQualityMVA)
                      SOA_COLUMN(double, pt),
                      SOA_COLUMN(double, eta),
                      SOA_COLUMN(float, phi),
                      SOA_COLUMN(float, dzBS),   // dz and dxy to the beam spot
                      SOA_COLUMN(float, dxyBS),
                      SOA_COLUMN(float, chi2),
                      SOA_COLUMN(float, ndof),
                      SOA_COLUMN(int16_t, npixBarrel),
                      SOA_COLUMN(int16_t, npixEndcap),
                      SOA_COLUMN(float, outerRadius),
                      SOA_COLUMN(float, outerZ),
                      // track parameters for dz to a vertex (TOF PID), as reco::TrackBase::dz
                      SOA_COLUMN(double, vx),
                      SOA_COLUMN(double, vy),
                      SOA_COLUMN(double, vz),
                      SOA_COLUMN(double, px),
                      SOA_COLUMN(double, py),
                      SOA_COLUMN(double, pz),
                      SOA_COLUMN(double, dzError))

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

  GENERATE_SOA_LAYOUT(MtdLayerLayout,
                      SOA_COLUMN(float, ref),        // BTL: mean module radius; ETL: mean module z
                      SOA_COLUMN(float, spread),     // max |r - ref| (BTL) or |z - ref| (ETL) of the modules
                      SOA_COLUMN(float, maxHitErr))  // max sqrt(hit local error) in the layer

  GENERATE_SOA_BLOCKS(MtdEventInputBlocks,
                      SOA_BLOCK(tracks, MtdTrackInLayout),
                      SOA_BLOCK(hits, MtdHitInLayout),
                      SOA_BLOCK(cells, MtdCellLayout),
                      SOA_BLOCK(layers, MtdLayerLayout),
                      SOA_BLOCK(vtxTracks, reco::TrackForVertexSoALayout))

  using MtdTrackInSoA = MtdTrackInLayout<>;
  using MtdHitInSoA = MtdHitInLayout<>;
  using MtdCellSoA = MtdCellLayout<>;
  using MtdLayerSoA = MtdLayerLayout<>;
  using MtdEventInputSoA = MtdEventInputBlocks<>;

  constexpr int kNLayers = 5;

}  // namespace mtdtiming

#endif
