#ifndef DataFormats_MtdTimingSoA_interface_MtdTimingHost_h
#define DataFormats_MtdTimingSoA_interface_MtdTimingHost_h

#include "DataFormats/MtdTimingSoA/interface/MtdTimingSoA.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"

namespace mtdtiming {
  using MtdTimingHost = PortableHostCollection<MtdTimingSoA>;
}

#endif
