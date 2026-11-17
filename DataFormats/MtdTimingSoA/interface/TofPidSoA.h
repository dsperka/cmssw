#ifndef DataFormats_MtdTimingSoA_interface_TofPidSoA_h
#define DataFormats_MtdTimingSoA_interface_TofPidSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace mtdtiming {

  // One row per input track (same order as the input track collection): the output of the TOF
  // particle identification (TOFPIDProducer), computed on the device. Same defaults as TOFPIDProducer.
  GENERATE_SOA_LAYOUT(TofPidLayout,
                      SOA_COLUMN(float, t0),
                      SOA_COLUMN(float, sigmat0),
                      SOA_COLUMN(float, t0safe),
                      SOA_COLUMN(float, sigmat0safe),
                      SOA_COLUMN(float, probPi),
                      SOA_COLUMN(float, probK),
                      SOA_COLUMN(float, probP))

  using TofPidSoA = TofPidLayout<>;

}  // namespace mtdtiming

#endif
