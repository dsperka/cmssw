#ifndef DataFormats_VertexSoA_interface_TrackForVertexHostCollection_h
#define DataFormats_VertexSoA_interface_TrackForVertexHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/VertexSoA/interface/TrackForVertexSoA.h"

// SoA with x, y, z, id fields in host memory
namespace reco {
  using TrackForVertexHostCollection = PortableHostCollection<TrackForVertexSoA>;
}  // namespace reco

#endif  // DataFormats_VertexSoA_interface_TrackForVertexHostCollection_h
