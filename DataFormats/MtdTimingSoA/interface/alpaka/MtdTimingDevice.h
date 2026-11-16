#ifndef DataFormats_MtdTimingSoA_interface_alpaka_MtdTimingDevice_h
#define DataFormats_MtdTimingSoA_interface_alpaka_MtdTimingDevice_h

#include "DataFormats/MtdTimingSoA/interface/MtdTimingHost.h"
#include "DataFormats/MtdTimingSoA/interface/MtdTimingSoA.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming {
  using namespace ::mtdtiming;
  using MtdTimingDevice = PortableCollection<MtdTimingSoA>;
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(mtdtiming::MtdTimingDevice, mtdtiming::MtdTimingHost);

#endif
