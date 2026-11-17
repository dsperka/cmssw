#ifndef DataFormats_VertexSoA_interface_alpaka_VertexTrackDeviceCollection_h
#define DataFormats_VertexSoA_interface_alpaka_VertexTrackDeviceCollection_h

#include "DataFormats/VertexSoA/interface/VertexTrackHostCollection.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/VertexSoA/interface/VertexTrackSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using ::reco::VertexTrackHostCollection;

  using VertexTrackDeviceCollection = PortableCollection<::reco::VertexTrackSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // DataFormats_VertexSoA_interface_alpaka_VertexTrackDeviceCollection_h
