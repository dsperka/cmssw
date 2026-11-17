#ifndef DataFormats_MtdTimingSoA_interface_MtdTimingSoA_h
#define DataFormats_MtdTimingSoA_interface_MtdTimingSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace mtdtiming {

  // One row per input track (same order as the input track collection).
  // Times in ns, lengths in cm. Defaults (no MTD time) follow the legacy TrackExtenderWithMTD:
  //   chi2 = -1, pathLength = -1, tmtd = 0, sigmatmtd = -1, tof = 0, sigmaTof = -1, t0 = 0, sigmat0 = -1,
  //   beta = 0, position = (0,0,0).
  GENERATE_SOA_LAYOUT(MtdTimingLayout,
                      SOA_COLUMN(float, btlChi2),       // position chi2 of the chosen BTL hit
                      SOA_COLUMN(float, btlTimeChi2),   // time chi2 of the chosen BTL hit
                      SOA_COLUMN(float, etlChi2),       // best of the chosen ETL hits
                      SOA_COLUMN(float, etlTimeChi2),
                      SOA_COLUMN(int32_t, btlHit),      // index of the chosen hit in the input MTDTrackingDetSetVector data, -1 if none
                      SOA_COLUMN(int32_t, etlHit1),     // inner ETL hit
                      SOA_COLUMN(int32_t, etlHit2),     // outer ETL hit
                      SOA_COLUMN(int8_t, timed),        // 1 if the track has an MTD time
                      SOA_COLUMN(float, pathLength),
                      SOA_COLUMN(float, tmtd),
                      SOA_COLUMN(float, sigmatmtd),
                      SOA_COLUMN(float, mtdPosX),
                      SOA_COLUMN(float, mtdPosY),
                      SOA_COLUMN(float, mtdPosZ),
                      SOA_COLUMN(float, tofPi),
                      SOA_COLUMN(float, tofK),
                      SOA_COLUMN(float, tofP),
                      SOA_COLUMN(float, sigmaTofPi),
                      SOA_COLUMN(float, sigmaTofK),
                      SOA_COLUMN(float, sigmaTofP),
                      SOA_COLUMN(float, t0),
                      SOA_COLUMN(float, sigmat0),
                      SOA_COLUMN(float, beta),
                      SOA_COLUMN(float, sigmaBeta),
                      SOA_COLUMN(float, mva))  // MTD track-quality BDT (mtdTrackQualityMVA:mtdQualMVA), -1 if not computed

  using MtdTimingSoA = MtdTimingLayout<>;
  using MtdTimingView = MtdTimingSoA::View;
  using MtdTimingConstView = MtdTimingSoA::ConstView;

}  // namespace mtdtiming

#endif
