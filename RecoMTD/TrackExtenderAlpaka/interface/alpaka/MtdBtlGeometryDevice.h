#ifndef RecoMTD_TrackExtenderAlpaka_interface_alpaka_MtdBtlGeometryDevice_h
#define RecoMTD_TrackExtenderAlpaka_interface_alpaka_MtdBtlGeometryDevice_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoMTD/TrackExtenderAlpaka/interface/MtdBtlGeometryHost.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming {
  using namespace ::mtdtiming;
  using MtdBtlGeometryDevice = PortableCollection<MtdBtlGeometrySoA>;
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming

#endif
