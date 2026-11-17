#ifndef DataFormats_MtdTimingSoA_interface_TofPidHost_h
#define DataFormats_MtdTimingSoA_interface_TofPidHost_h

#include "DataFormats/MtdTimingSoA/interface/TofPidSoA.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"

namespace mtdtiming {
  using TofPidHost = PortableHostCollection<TofPidSoA>;
}

#endif
