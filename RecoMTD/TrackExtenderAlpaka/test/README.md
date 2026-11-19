# RecoMTD/TrackExtenderAlpaka tests

How to reproduce the checks of the Alpaka MTD chain: the MTD track timing (no refit), and on top of it the
3D vertices (DA in blocks, both fits from one clustering, vertex time) and the TOF PID on the device.

| directory | what it shows |
|---|---|
| `devicechain/` | the device chain gives exactly the same numbers as the host code it replaces (100-event skim) |
| `validation/` | physics effect on the 4D vertices and the MTD track validation, legacy vs Alpaka (1000 events) |
| `timing/` | time per event of the full reconstruction and per module, legacy vs Alpaka (300 events, 1 thread) |
| `refit/` | what the legacy refit changes: momentum, path length and t0 against the truth (100-event skim) |
| `testMtdTimingMath.cpp`, `MtdTimingCompare.cc`, `compare_cfg.py` | unit test of the track math; per-track legacy vs Alpaka dump |

## Setup

```
cmsrel CMSSW_20_1_0_pre5 && cd CMSSW_20_1_0_pre5/src && cmsenv
# add the branch, then
scram b -j 8      # -j 20 can run out of memory on the Alpaka code
```

Input: RelValTTbar_14TeV, D128, PU200, GEN-SIM-RECO of CMSSW_20_1_0_pre3, listed in
`validation/relvalTTbar_PU200_D128_20_1_0_pre3.txt` (10 files of 100 events). The scripts read them with
xrootd by default; with local copies, pass `file:/path/to/` as the file prefix. These are 17 GB files:
copying them is worth it, and needed for the timing (the page cache is dropped between jobs).
The jobs run on CPU (`accelerators = ["cpu"]`).

The RelVal digis are made with the 1000 fb⁻¹ detector aging, so the reconstruction must apply it too, as the
RelVal step3 does (`SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000`). `mtdValidation_cfg.py`
does this by default (option `aging`, 1000; `aging=0` for none), and the skim, timing and refit configurations
inherit it. Without it the BTL hits get the unaged light output and time walk: hit times about 0.7 ns late
(more at low energy) and energies 1.4× too high, which splits 4D vertices.

## The customisation

`python/customizeMtdTimingAlpaka.py`:

- `customizeMtdTimingAlpaka(process)`: the Alpaka MTD timing in place of `trackExtenderWithMTD` (same label and
  product names), and the quality BDT on the device (`mtdTrackQualityMVA` becomes an alias).
- `customizeMtdDeviceVertices(process)` (after the first one): `unsortedOfflinePrimaryVertices` (no beam
  spot and `WithBS`) and `tofPID4DnoPID` on the device, from the same packed input. Same labels and product
  names, so the sorters, the 4D vertexing and the validation run unchanged.

In `validation/mtdValidation_cfg.py` and `timing/mtdTiming_cfg.py` they are the options `timing=alpaka` and
`vertices=alpaka`.

## How the Alpaka code works, and how it differs from legacy

### Data flow

```
host                                  device                                        host
generalTracks, MTD rechits,  ──1 copy──>  MTD timing + BDT ─┬─> 3D vertices ──> vertex time ──┬─> converters ─> same products
beam spot, (BTL geometry:                                   └────────────────> TOF PID ───────┘   as legacy
once per IOV; BDT: once per job)
```

- `MtdEventInputProducer` (host) packs everything the device needs from the event into one product: each
  track at its outermost hit (position, momentum, 5x5 error) plus the numbers the BDT and TOF PID use, the
  MTD hits sorted into search cells, a few per-layer numbers, and the tracks selected for vertexing, sorted
  in z (same selection as the standard vertexing). The framework copies it to the device once; on CPU there
  is no copy. Small numbers (cuts, beam spot x/y/width) go to the kernels as arguments.
- The BTL geometry (rods and modules as planes) is an event setup product (`MtdBtlGeometryESProducer`),
  copied once per geometry. The BDT trees are flattened into tables and copied once per job.
- The device results come back as SoA tables and are turned into the legacy products by host converters
  (`MtdTimingToValueMaps`, `RecoVertexFromSoAProducer`, `TofPidToValueMaps`), under the legacy labels.
  Today that is seven device-to-host copies.

### MTD track timing (`MtdTimingProducer`, replaces `TrackExtenderWithMTD`)

Three kernels:

1. **Hit matching**, one thread per track and MTD layer (BTL, two ETL discs per side). The track is moved
   from its outermost tracker hit to the layer along a helix in the field at its start point (as the legacy
   analytical propagator). The candidate hits are found with the legacy navigation rules: in BTL the closest
   rod in phi (and the next one when the track is near the edge), then the modules along z; in ETL the modules
   around the crossing point. Every candidate gets the legacy position chi2 (track and hit errors) and time
   chi2, the best one is kept by the legacy ordering, and only the best one is tested against the cuts.
2. **Timing**, one thread per track: combine the matched layers as legacy `buildTrack` does (one BTL hit, or
   one or two ETL hits; BTL+ETL is rejected as legacy does), path length, time of flight for pion, kaon and
   proton, t0 and its error, beta.
3. **Quality BDT**, one thread per track: the `MTDTrackQualityMVA` BDT with the same inputs, cuts and trees.

The one real difference: **there is no refit.** Legacy refits each track with its MTD hit (Kalman fit) and
uses the refitted momentum and trajectory for the path length and the time of flight. The Alpaka version
keeps the original track: the path in the tracker is computed from the helix between the point of closest
approach and the outermost hit, and the momentum is the original one. Consequences:

- the output tracks are copies of the input tracks with the MTD hit added (`trackClones`), not refitted;
- σ(TOF) is 1–15% larger track by track (median; the refit makes the quoted momentum error 1–15% smaller);
- the path length differs by more than 1 cm for 27% of the BTL tracks and 14% of the ETL tracks (mostly at
  low momentum), so the t0 values differ, which in turn moves the 4D vertexing a little (+0.6 fake and +0.6
  split vertices per event out of ~130, real vertices unchanged).

Against the truth (section 4), the refit does not make the momentum more accurate: it is about the same in
BTL, worse in ETL, and the tails get larger in both. Where the path lengths differ by more than 5 cm (8% of
the BTL tracks, mostly below 1 GeV), the Alpaka t0 is the one closer to the true time for 74% of the pions. So no
refit (or an approximate one) is planned.
On CPU the MTD time, path and hit choice were checked against legacy track by track (`MtdTimingCompare.cc`);
the BDT output is identical to `MTDTrackQualityMVAProducer` on the same inputs.

### 3D vertices (`PrimaryVertexProducerPortable`, replaces `unsortedOfflinePrimaryVertices`)

| | legacy | Alpaka |
|---|---|---|
| clustering | `DA_vect`: deterministic annealing over all tracks of the event at once | DA in blocks: tracks sorted in z, cut into blocks of 512 with 50% overlap, annealing in each block in parallel, then the vertices of overlapping blocks are merged and each track is assigned to one vertex |
| fit | `AdaptiveVertexFitter`, once per collection | weighted mean (as `WeightedMeanFitter`), iterated with track weights |
| no BS and WithBS | two separate clusterings and fits | one clustering, two fits |
| vertex time | `fromTracksPID` | the same algorithm on the device, same numbers |

So the device vertices are not expected to be identical to legacy; the effect is in the validation table
below (about +0.2 real, +0.3 fake, +0.3 split vertices per event compared to host vertices).

### TOF PID (`TofPidProducer`, replaces `tofPID4DnoPID`)

A direct port of `TOFPIDProducer`, one thread per track: find the track's vertex (from the device 3D
vertices), compare the track time with the vertex time under the pion, kaon and proton hypotheses, return
t0, σ(t0) and the three probabilities. Same parameters, same seven output maps, identical values on CPU.

### Exactness on CPU

Where the Alpaka code ports a legacy algorithm (BDT, vertex time, TOF PID) it gives the same numbers on the
CPU backend, bit for bit (`devicechain/`). To get there the device code uses the same float/double steps as
the host code (for example the same fast exp on the host and the C `sqrt(double)` in the BDT). On a GPU,
`exp` and `sqrt` come from the GPU math library, so last-bit differences are expected there; the GPU has not
been tested yet.

## 1. Exact comparison with the host code (`devicechain/`)

Make the skim once (full reconstruction of 100 events, keeps tracks, MTD clusters and rechits, beam spot):

```
cmsRun devicechain/skim_cfg.py inputFiles=file:/path/to/01313427-f364-4b67-9190-ea6bafe583b6.root \
    maxEvents=100 output=skim.root threads=8
```

Then run the chain (`chain_cfg.py`, MTD rechits rebuilt from the clusters) and compare:

```
cd devicechain
cmsRun chain_cfg.py inputFiles=file:skim.root mode=dev2 output=dev2.root threads=8
python3 compare_pairs.py dev2.root
```

`mode=dev2` runs the device chain (packed input, MTD timing with the BDT, one clustering with both fits and
the vertex time, TOF PID) and, next to it, the host `MTDTrackQualityMVAProducer` and `TOFPIDProducer` on the
same vertices. `compare_pairs.py` compares the BDT output and the 7 TOF PID maps value by value.
Expected: `TOTAL differing values: 0`.

Other modes: `new` (packed input, two separate vertex producers) and `tracks` (the vertices read the separate
`PortableTrackSoAProducer` tables). `compare_chain.py a.root b.root` compares two outputs leaf by leaf
(every value map and vertex): `new` and `tracks` must agree exactly, and a rerun must agree with itself.

## 2. Physics validation (`validation/`)

```
mkdir work && cd work
../validation/run_validation.sh 1000 8 ../validation/relvalTTbar_PU200_D128_20_1_0_pre3.txt file:/path/to
```

Runs the reconstruction twice in parallel (legacy, then `timing=alpaka`) with `MtdTracksValidation` and
`Primary4DVertexValidation`, harvests, writes `summary.html/json` (vertices per event, primary-vertex
track composition, vertex time resolution) and the plot galleries `mtd_{Vertices,Tracks}_gallery.html`.
For the device vertices and TOF PID: `NEW_OPTS="vertices=alpaka" NEW_NAME="..." run_validation.sh ...`.

Results on 1000 events (the 10 files), 8 threads. Columns: legacy TrackExtenderWithMTD (with refit) and host
vertices; Alpaka MTD timing (no refit) and host vertices; Alpaka MTD timing with the device 3D vertices and
TOF PID (`vertices=alpaka`). 4D vertices (`offlinePrimaryVertices4D`), from `mtd_summary.py`.

| per event | legacy | Alpaka timing | + device vertices, TOF PID |
|---|---|---|---|
| simulated vertices | 189.6 | 189.6 | 189.6 |
| reconstructed, all | 142.0 | 142.5 | 142.9 |
| reconstructed, pass ndof | 119.8 | 120.2 | 120.7 |
| real | 108.99 | 109.15 | 109.37 |
| fake | 10.81 | 11.05 | 11.32 |
| split | 5.24 | 5.45 | 5.70 |
| other fake | 5.57 | 5.59 | 5.62 |
| real / (real + fake) | 0.910 | 0.908 | 0.906 |
| real / simulated | 0.575 | 0.576 | 0.577 |

| leading vertex, mean fraction of its tracks (lower is better) | legacy | Alpaka timing | + device vertices, TOF PID |
|---|---|---|---|
| pileup tracks, by multiplicity | 0.2302 | 0.2280 | 0.2297 |
| pileup tracks, sum pT | 0.0878 | 0.0849 | 0.0859 |
| pileup tracks, weighted | 0.0193 | 0.0178 | 0.0185 |
| fake tracks, by multiplicity | 0.0484 | 0.0479 | 0.0474 |
| secondary tracks, weighted | 0.0619 | 0.0626 | 0.0621 |

| vertex time | legacy | Alpaka timing | + device vertices, TOF PID |
|---|---|---|---|
| RMS of t(rec) − t(sim), all vertices [ps] | 49.7 | 50.5 | 50.3 |
| RMS, signal vertex [ps] | 12.8 | 12.7 | 13.4 |
| fraction within 30 ps | 0.651 | 0.644 | 0.643 |
| vertices with a time (1000 events) | 108986 | 109154 | 109374 |

Without the aging customisation (earlier runs) all three columns had about 14.5 split and 20.4 fake vertices
per event, a signal-vertex time RMS of 20 ps and 47% of the vertices within 30 ps; the Alpaka timing then
gave 0.6 more fake and split vertices per event than legacy, against about 0.2 now.

The device 3D vertices are DA in blocks with a weighted-mean fit, not the host clusterizer and fitter; that
is where the small differences in the last column come from. Two legacy runs give identical numbers.

## 3. Timing (`timing/`)

```
mkdir work && cd work
../timing/run_timing.sh 300 1 ../validation/relvalTTbar_PU200_D128_20_1_0_pre3.txt file:/path/to 2
```

- `mtdTiming_cfg.py`: the reconstruction of `mtdValidation_cfg.py` without validation, DQM or output file.
  A `GenericConsumer` asks for the RECOSIM products made in the job, so the same modules run as in a RECO
  step and nothing is written. FastTimerService (per module, `<job>.json`) and ThroughputService (events/s,
  first 5 events skipped; `MTD_TIMING_SKIP` changes it). Threads = streams.
- `readOnly_cfg.py`: reads only the input branches the reconstruction reads (from the job report of a short
  job, `branches.txt`) and nothing else: the most the disk can deliver.
- `run_timing.sh`: read-only job, then legacy, Alpaka timing, Alpaka timing + device vertices and TOF PID,
  one after the other on one pinned core. Before each job: input files dropped from the page cache
  (`drop_page_cache.py`, no root needed), wait for an idle machine. `monitor.sh` logs load, core frequency,
  temperature, memory and disk rate every 15 s. At the end `timing_summary.py` prints the tables.

Results (measured before the aging customisation was added; to be redone), 300 events, 1 thread, Intel Xeon W-2155 (turbo on, governor powersave), input on a local HDD:

| throughput (events 6–300) | events/s | s/event |
|---|---|---|
| read only | 0.602 | 1.66 |
| legacy | 0.02822 | 35.44 |
| Alpaka timing | 0.02909 | 34.38 (−3.0%) |
| + device vertices, TOF PID | 0.02947 | 33.94 (−4.2%) |

| ms per event (thread time) | legacy | Alpaka timing | + device vertices, TOF PID |
|---|---|---|---|
| MTD timing + quality BDT | 1128 | 204 | 202 |
| 3D vertices (both fits) | 710 | 706 | 291 |
| TOF PID (tofPID4DnoPID) | 33 | 34 | 2.7 |
| total | 1872 | 943 | 495 |
| all other modules | 31708 | 31588 | 31594 |

The input is not a limit (reading takes 1.66 s of 34–35 s). The unchanged modules agree within 0.4%, which is
the run-to-run noise of this machine (the core clock moves with turbo and temperature; mean 3790–3854 MHz).

## 4. Refit study (`refit/`)

Legacy refits each track with its MTD hit; the Alpaka version does not. This study checks what the refit
changes, track by track, and which version is closer to the truth. Input: the skim of section 1.

```
mkdir work && cd work
../refit/run_refit_study.sh /path/to/skim.root 8
```

- the Alpaka MTD timing on the skim (`devicechain/chain_cfg.py mode=new`) → `chain.root`;
- `make_table.py` (FWLite): one row per track with an MTD time in both legacy and Alpaka: original and
  refitted momentum and error at the point of closest approach (the refitted tracks have the same index as
  `generalTracks`), the legacy and Alpaka value maps, and the closest simulated particle (same charge,
  ΔR < 0.005, |Δpt/pt| < 0.3) with its momentum, type and production time → `table.npz`;
- `analyze_refit.py` → `refit.txt`: what the refit does to the momentum, refit vs truth, legacy vs Alpaka
  inputs of the time of flight;
- `analyze_truth.py` → `truth.txt`: t0 (pion mass) minus the true production time, for simulated pions from
  in-time collisions produced at the beam line.

Results, 100 events: 314k tracks with an MTD time in both, 126k matched to a simulated particle, 92k of them
prompt in-time pions.

**The refit barely changes the quoted momentum error, and legacy σ(TOF) is exactly the refitted one.** The
refitted σ(p) is 0.94–0.99 of the original in BTL and 0.85–0.96 in ETL (median per track); the σ(p) implied
by legacy σ(TOF) equals the refitted σ(p) for every track, and the legacy momentum map is the refitted
momentum. So Alpaka σ(TOF) / legacy (median per track) is 1.01–1.06 in BTL and 1.02–1.15 in ETL.

**The refitted momentum is not closer to the truth.** |p/p_true − 1| in %, original → refit:

| p [GeV], region | 68% of tracks within | 95% within | tracks off by more than 10% |
|---|---|---|---|
| 0.7–1, BTL | 1.05 → 1.09 | 3.7 → 7.9 | 2.9 → 4.5 |
| 1–1.5, BTL | 0.99 → 1.09 | 3.2 → 5.6 | 2.0 → 3.1 |
| 2–3, BTL | 0.93 → 0.95 | 2.6 → 3.7 | 0.4 → 0.9 |
| 1–1.5, ETL | 2.3 → 3.8 | 21 → 26 | 15 → 21 |
| 2–3, ETL | 2.8 → 4.0 | 19 → 36 | 10 → 21 |
| above 5, ETL | 3.1 → 3.6 | 8.7 → 20 | 3.5 → 12 |

The pulls ((p − p_true)/σ(p)) also get wider after the refit. Matching is by direction only, so part of the
ETL tails are wrong matches; they are the same tracks for the original and the refit. The true momentum is
taken at production, before energy loss. `make_plots.py` draws these (`plots/`).

**Path length.** The medians agree, but the path differs by more than 1 cm (33 ps) for 25% of the BTL and
14% of the ETL tracks, and by more than 5 cm for 7% and 3% (18% of BTL tracks at 0.7–1 GeV). Legacy and
Alpaka pick the same MTD hit for 97.8% of the pions.

**t0 against the true time** (pions, same MTD hit):

| path difference | tracks | median \|t0 − t_true\| legacy | Alpaka | Alpaka closer |
|---|---|---|---|---|
| BTL, 1–5 cm | 5888 | 110 ps | 123 ps | 43% |
| BTL, more than 5 cm (median p 0.86 GeV) | 2511 | 950 ps | 131 ps | 86% |
| ETL, 1–5 cm | 1837 | 268 ps | 266 ps | 52% |
| ETL, more than 5 cm | 455 | 428 ps | 396 ps | 52% |

Width of t0 − t_true (half of the central 68%), legacy / Alpaka:

| p [GeV] | 0.7–1 | 1–1.5 | 1.5–2 | 2–3 | 3–5 | above 5 |
|---|---|---|---|---|---|---|
| BTL [ps] | 374 / 85 | 81 / 83 | 80 / 84 | 71 / 74 | 65 / 66 | 56 / 57 |
| ETL [ps] | 60 / 60 | 52 / 53 | 60 / 59 | 92 / 100 | 66 / 68 | 49 / 48 |

Where the two path lengths disagree a lot (mostly BTL below 1 GeV), it is the legacy refitted trajectory that
is wrong. For small differences (1–5 cm) legacy is slightly closer. The medians of t0 − t_true are +1 ps in BTL
and +4 ps in ETL. σ(t0) is large in both (median 1.23 ns in BTL, 0.29 ns in ETL, against 41 and 26 ps for the
MTD hit), so it must include the unknown-mass term.
