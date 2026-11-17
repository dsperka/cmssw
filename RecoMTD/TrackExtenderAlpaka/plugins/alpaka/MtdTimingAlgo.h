#ifndef RecoMTD_TrackExtenderAlpaka_plugins_alpaka_MtdTimingAlgo_h
#define RecoMTD_TrackExtenderAlpaka_plugins_alpaka_MtdTimingAlgo_h

#include "DataFormats/MtdTimingSoA/interface/alpaka/MtdTimingDevice.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoMTD/TrackExtenderAlpaka/interface/MtdQualityBdt.h"
#include "RecoMTD/TrackExtenderAlpaka/interface/MtdTimingInputSoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming {

  using RecordDevice = PortableCollection<::mtdtiming::MtdRecordSoA>;

  // Three kernels:
  //   matchHits    one thread per (track, layer): best hit in that layer
  //   finishTiming one thread per track: combine the layers, path length, TOF, t0
  //   quality      one thread per track: the MTD track-quality BDT
  void runMtdTiming(Queue& queue,
                    ::mtdtiming::MatchParams const& params,
                    ::mtdtiming::MtdEventInputSoA::ConstView input,
                    ::mtdtiming::MtdBtlGeometrySoA::ConstView geometry,
                    ::mtdtiming::MtdBdtSoA::ConstView bdt,
                    RecordDevice& records,
                    MtdTimingDevice& out);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::mtdtiming

#endif
