#ifndef DataFormats_VertexSoA_interface_VertexTrackHostCollection_h
#define DataFormats_VertexSoA_interface_VertexTrackHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/VertexSoA/interface/VertexTrackSoA.h"

namespace reco {
  using VertexTrackHostCollection = PortableHostCollection<VertexTrackSoA>;
}  // namespace reco

#endif  // DataFormats_VertexSoA_interface_VertexTrackHostCollection_h
