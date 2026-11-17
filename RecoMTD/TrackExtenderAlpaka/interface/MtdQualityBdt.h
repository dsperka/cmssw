#ifndef RecoMTD_TrackExtenderAlpaka_interface_MtdQualityBdt_h
#define RecoMTD_TrackExtenderAlpaka_interface_MtdQualityBdt_h

// The MTD track-quality BDT (RecoMTD/TimingIDTools MTDTrackQualityMVA) as flat tables, so that it can
// be evaluated on the device. Same trees and same walk as GBRForest/GBRTree.
//   trees:     one row per tree, its first node and first leaf response
//   nodes:     cut variable, cut value, left and right child (a child <= 0 is the leaf -child)
//   responses: leaf responses

#include <cmath>
#include <cstdint>

#include <alpaka/core/Common.hpp>

#include "DataFormats/SoATemplate/interface/SoABlocks.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace mtdtiming {

  GENERATE_SOA_LAYOUT(MtdBdtTreeLayout,
                      SOA_COLUMN(int32_t, firstNode),
                      SOA_COLUMN(int32_t, firstResponse),
                      SOA_SCALAR(double, initialResponse))
  GENERATE_SOA_LAYOUT(MtdBdtNodeLayout,
                      SOA_COLUMN(uint8_t, cutIndex),
                      SOA_COLUMN(float, cutVal),
                      SOA_COLUMN(int32_t, left),
                      SOA_COLUMN(int32_t, right))
  GENERATE_SOA_LAYOUT(MtdBdtResponseLayout, SOA_COLUMN(float, response))

  // sizes: trees, nodes, responses
  GENERATE_SOA_BLOCKS(MtdBdtBlocks,
                      SOA_BLOCK(trees, MtdBdtTreeLayout),
                      SOA_BLOCK(nodes, MtdBdtNodeLayout),
                      SOA_BLOCK(responses, MtdBdtResponseLayout))
  using MtdBdtSoA = MtdBdtBlocks<>;

  // the 17 inputs, in the order of MTDTRACKQUALITYMVA_VARS
  constexpr int kNBdtVars = 17;

  // GBRForest::GetResponse
  ALPAKA_FN_HOST_ACC inline double bdtResponse(MtdBdtSoA::ConstView bdt, float const* vars) {
    auto const trees = bdt.trees();
    auto const nodes = bdt.nodes();
    auto const responses = bdt.responses();
    double response = trees.initialResponse();
    for (int32_t t = 0; t < trees.metadata().size(); ++t) {
      int32_t const first = trees[t].firstNode();
      int32_t index = 0;
      do {
        auto const n = nodes[first + index];
        index = vars[n.cutIndex()] > n.cutVal() ? n.right() : n.left();
      } while (index > 0);
      response += responses[trees[t].firstResponse() - index].response();
    }
    return response;
  }

  // MTDTrackQualityMVA: GBRForest::GetGradBoostClassifier, returned as float by TMVAEvaluator, then
  // mapped to a probability in [0, 1]
  ALPAKA_FN_HOST_ACC inline float bdtQuality(MtdBdtSoA::ConstView bdt, float const* vars) {
    double const response = bdtResponse(bdt, vars);
    float const value = 2.0 / (1.0 + std::exp(-2.0 * response)) - 1;
    // the host code calls the C sqrt(double) here
    return 1. / (1 + std::sqrt(double(2 / (1 + value) - 1)));
  }

}  // namespace mtdtiming

#endif
