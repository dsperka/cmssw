#ifndef DataFormats_MtdTimingSoA_interface_MtdEventInputHost_h
#define DataFormats_MtdTimingSoA_interface_MtdEventInputHost_h

#include "DataFormats/MtdTimingSoA/interface/MtdEventInputSoA.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"

namespace mtdtiming {
  // sizes: tracks, hits, cells, layers, vtxTracks
  using MtdEventInputHost = PortableHostCollection<MtdEventInputSoA>;
}

#endif
