#ifndef DataFormats_VertexSoA_interface_VertexSoA_h
#define DataFormats_VertexSoA_interface_VertexSoA_h

#include <Eigen/Core>
#include <Eigen/Dense>
#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace reco {
  GENERATE_SOA_LAYOUT(
      VertexSoALayout,
      SOA_COLUMN(float, x),
      SOA_COLUMN(float, y),
      SOA_COLUMN(float, z),
      SOA_COLUMN(float, t),

      SOA_COLUMN(float, errx),
      SOA_COLUMN(float, erry),
      SOA_COLUMN(float, errz),
      SOA_COLUMN(float, errt),

      SOA_COLUMN(float, chi2),
      SOA_COLUMN(float, ndof),
      SOA_COLUMN(int, ntracks),

      // The per-vertex scratch of the clusterizer (rho, aux1, aux2, sw, swE) is not part of this layout: it lives in a
      // device-only collection, see DAInBlocksScratch.h in the plugin

      // The tracks of the vertex are the entries first_track to first_track+ntracks-1 of the VertexTrackSoA that goes
      // with this collection
      SOA_COLUMN(int, first_track),

      SOA_COLUMN(bool, isGood),
      SOA_COLUMN(int, order),

      // Use entries for blocks
      // When running in paralell, i.e. DA in blocks, the entries in the column can be used to save the number of vertex on each block during clustering
      SOA_COLUMN(int32_t, nV))

  using VertexSoA = VertexSoALayout<>;
}  // namespace reco

#endif  // DataFormats_VertexSoA_interface_VertexSoA_h
