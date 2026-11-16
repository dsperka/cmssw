#ifndef RecoMTD_TrackExtenderAlpaka_plugins_alpaka_MtdTimingAlgo_h
#define RecoMTD_TrackExtenderAlpaka_plugins_alpaka_MtdTimingAlgo_h

#include "DataFormats/MtdTimingSoA/interface/alpaka/MtdTimingDevice.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoMTD/TrackExtenderAlpaka/interface/MtdTimingInputSoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming {

  using TrackInDevice = PortableCollection<::mtdtiming::MtdTrackInSoA>;
  using HitInDevice = PortableCollection<::mtdtiming::MtdHitInSoA>;
  using CellDevice = PortableCollection<::mtdtiming::MtdCellSoA>;
  using RecordDevice = PortableCollection<::mtdtiming::MtdRecordSoA>;
  using PlaneDevice = PortableCollection<::mtdtiming::MtdPlaneSoA>;

  // Two kernels:
  //   matchHits    one thread per (track, layer): best hit in that layer
  //   finishTiming one thread per track: combine the layers, path length, TOF, t0
  void runMtdTiming(Queue& queue,
                    ::mtdtiming::MatchParams const& params,
                    TrackInDevice const& tracks,
                    HitInDevice const& hits,
                    CellDevice const& cells,
                    PlaneDevice const& rods,
                    PlaneDevice const& modules,
                    RecordDevice& records,
                    MtdTimingDevice& out);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming

#endif
