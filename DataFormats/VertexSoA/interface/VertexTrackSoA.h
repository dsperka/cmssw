#ifndef DataFormats_VertexSoA_interface_VertexTrackSoA_h
#define DataFormats_VertexSoA_interface_VertexTrackSoA_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace reco {
  // Tracks of the vertices of a VertexSoA, as one list: the tracks of vertex v are the entries first_track(v) to
  // first_track(v)+ntracks(v)-1 of this collection
  GENERATE_SOA_LAYOUT(VertexTrackSoALayout,
                      SOA_COLUMN(int, track_id),  // index in the TrackForVertexSoA product of the same producer
                      SOA_COLUMN(float, track_weight))

  using VertexTrackSoA = VertexTrackSoALayout<>;
}  // namespace reco

#endif  // DataFormats_VertexSoA_interface_VertexTrackSoA_h
