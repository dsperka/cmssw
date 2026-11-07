# PR 50766 review campaign: validation record

Starting point: branch `from-CMSSW_16_1_0_pre4-pv-fixes-on-latest-pr` at 70adb1ade92, which was validated on Alpaka CPU and GPU.
The review comments of cms-sw/cmssw PR 50766 (94 inline threads plus the conversation items) were applied
one validated batch at a time. Every batch was rebuilt and run on Alpaka CPU and GPU (91 events of
RelValZMM_14 PU150, `Zmumu_16X_local.txt`), harvested, and compared bin by bin with the baseline.

## Acceptance rule
- Compared: all 311 1D histograms of `PrimaryVertexV`, Alpaka CPU and GPU separately.
- EXACT: no histogram differs. TOLERANT: every one of six key metrics (efficiency, fake rate, merged, duplicate,
  gen-duplicate, PV purity) moved by at most 2e-3. Anything else: batch reverted.
- The 2e-3 threshold is a judgement call, not derived from statistics.
- Each batch is compared with the original baseline, so small differences introduced by the float conversion
  (b03a/b03b: 46 histograms on CPU, 152 on GPU) are carried forward unchanged in later batches.

## Result
- 18 batches accepted, 1 reverted (b06, `Acc1D`: CPU PV purity moved by +0.0105, over the 2e-3 threshold).
- b17 (implicit host-to-device copy) failed to build twice (my errors in how the token was declared), then was accepted.
- Final tree vs baseline: efficiency 0.7435 (CPU) / 0.7442 (GPU), fake rate 0.1500 / 0.1493, merged 0.1537 / 0.1543.
  All identical to the baseline to four digits.
- Dispositions: 74 accepted, 7 already fixed, 5 deferred, 4 not applicable, 4 rejected. See `REVIEW_LOG.md`.

## Timing (ms per event, vertex producer module, single run each, 91 events)
| Mode | pre-campaign | final tree |
|---|---|---|
| Legacy D.A. | 369 | 366 |
| D.A. In Blocks (C++) | 87 | 87 |
| D.A. In Blocks (Alpaka CPU) | 182 | 192 |
| D.A. In Blocks (Alpaka GPU) | 36 | 22 |

Single runs carry noise. The GPU gain is probably from the implicit host-to-device copy (b17) and the removed sync (b08),
but this was not isolated. GPU includes first-event warm-up.

## Contents
- `REVIEW_LOG.md`: accept/reject for every review comment and every batch.
- `verdicts.txt`, `compare/*.json`: per-batch comparison results (including repeated attempts).
- `patches/*.patch`: the diff of each batch.
- `galleries/`: full 345-plot galleries, same binning, pre-campaign code and final tree (self-contained HTML).
- `phone/`: the phone-sized pages that were shared (review log, baseline plots, final plots, with old and new binning).
- `merge_threshold_experiment_patch.py` and the section below.

## Histogram binning used for the rebinned galleries (not part of this commit's code)
The analyzer in `Validation/RecoVertex` booked the vertex-count histograms over 0-130, but this PU150 sample has about
150-230 generated vertices per event, so all `vs_PU` and `vs_NumVertices` plots were empty (everything in overflow).
Tested in a separate work area: `nPUmin = 100`, `nPUbins = 232` (bins of 2), and track-count bins starting at 1.
These edits live in a copy of `Validation/RecoVertex`; they are not in this commit.

## Merge-threshold experiment (not part of this commit's code)
Idea: the cross-block merge floor of `2 * vertexSize = 0.012 cm` over-merges at high track multiplicity, where the z resolution
is much smaller. Tested: merge if closer than `min(3 * combined sigma, 4 * vertexSize = 0.024 cm)` (C++ clusterizer and Alpaka arbitrator).
Results on the same 91 events (before -> after, final binning):

| | C++ In Blocks | Alpaka CPU | Alpaka GPU |
|---|---|---|---|
| efficiency | 0.741 -> 0.750 | 0.744 -> 0.750 | 0.744 -> 0.750 |
| fake rate | 0.130 -> 0.128 | 0.150 -> 0.149 | 0.149 -> 0.150 |
| merged | 0.153 -> 0.156 | 0.154 -> 0.156 | 0.154 -> 0.157 |
| duplicate | 0.026 -> 0.039 | 0.038 -> 0.049 | 0.040 -> 0.052 |
| reco vertices / event | +2.2 | +1.7 | +1.8 |

Legacy D.A. was unchanged (control). The new threshold gives more vertices and a little more efficiency, but raises the
duplicate rate by about 1.2 points and does not lower the merged rate, so it was not adopted. A significance sweep (4 and 5 sigma)
was proposed and has not been run.

## Open items
- b06 (`Acc1D`) not bisected; comments 25, 28, 45 remain rejected.
- Deferred design items: threads 5, 9, 12 (transient part of the SoA, common track/vertex data structures, per-vertex 1024-element vectors), conversation C58 (workflow exercising the code).
- Stale untracked directory `DataFormats/OfflineVertexSoA/` (old package name) is not part of this commit and should be deleted.
- Not pushed.
