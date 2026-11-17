#ifndef DataFormats_MtdTimingSoA_interface_alpaka_MtdEventInputDevice_h
#define DataFormats_MtdTimingSoA_interface_alpaka_MtdEventInputDevice_h

#include "DataFormats/MtdTimingSoA/interface/MtdEventInputHost.h"
#include "DataFormats/MtdTimingSoA/interface/MtdEventInputSoA.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming {
  using namespace ::mtdtiming;
  using MtdEventInputDevice = PortableCollection<MtdEventInputSoA>;
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(mtdtiming::MtdEventInputDevice, mtdtiming::MtdEventInputHost);

#endif
