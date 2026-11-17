#ifndef RecoMTD_TrackExtenderAlpaka_interface_MtdBtlGeometryHost_h
#define RecoMTD_TrackExtenderAlpaka_interface_MtdBtlGeometryHost_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "RecoMTD/TrackExtenderAlpaka/interface/MtdTimingInputSoA.h"

namespace mtdtiming {
  // event setup product; sizes: rods, modules, lookup (= modules)
  using MtdBtlGeometryHost = PortableHostCollection<MtdBtlGeometrySoA>;
}

#endif
