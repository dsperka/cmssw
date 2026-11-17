#ifndef DataFormats_MtdTimingSoA_interface_alpaka_TofPidDevice_h
#define DataFormats_MtdTimingSoA_interface_alpaka_TofPidDevice_h

#include "DataFormats/MtdTimingSoA/interface/TofPidHost.h"
#include "DataFormats/MtdTimingSoA/interface/TofPidSoA.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming {
  using namespace ::mtdtiming;
  using TofPidDevice = PortableCollection<TofPidSoA>;
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(mtdtiming::TofPidDevice, mtdtiming::TofPidHost);

#endif
