# RecoMTD/TrackExtenderAlpaka tests

How to reproduce the checks of the Alpaka MTD chain: the MTD track timing (no refit), and on top of it the
3D vertices (DA in blocks, both fits from one clustering, vertex time) and the TOF PID on the device.

| directory | what it shows |
|---|---|
| `devicechain/` | the device chain gives exactly the same numbers as the host code it replaces (100-event skim) |
| `validation/` | physics effect on the 4D vertices and the MTD track validation, legacy vs Alpaka (1000 events) |
| `timing/` | time per event of the full reconstruction and per module, legacy vs Alpaka (300 events, 1 thread) |
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
- σ(TOF) is about 35% larger (the momentum error of the original track is larger than after the refit);
- the t0 values differ slightly, which in turn moves the 4D vertexing a little (+0.6 fake and +0.6 split
  vertices per event out of ~130, real vertices unchanged).

A refit option for low-pt tracks (where it matters most) is not implemented.
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
| reconstructed, all | 158.2 | 158.8 | 159.7 |
| reconstructed, pass ndof | 130.4 | 131.0 | 131.6 |
| real | 109.96 | 110.03 | 110.25 |
| fake | 20.41 | 21.00 | 21.34 |
| split | 14.53 | 15.09 | 15.42 |
| other fake | 5.87 | 5.90 | 5.92 |
| real / (real + fake) | 0.843 | 0.840 | 0.838 |
| real / simulated | 0.580 | 0.580 | 0.581 |

| leading vertex, mean fraction of its tracks (lower is better) | legacy | Alpaka timing | + device vertices, TOF PID |
|---|---|---|---|
| pileup tracks, by multiplicity | 0.2114 | 0.2141 | 0.2131 |
| pileup tracks, sum pT | 0.0796 | 0.0806 | 0.0786 |
| pileup tracks, weighted | 0.0278 | 0.0285 | 0.0280 |
| fake tracks, by multiplicity | 0.0463 | 0.0456 | 0.0449 |
| secondary tracks, weighted | 0.0570 | 0.0578 | 0.0569 |

| vertex time | legacy | Alpaka timing | + device vertices, TOF PID |
|---|---|---|---|
| RMS of t(rec) − t(sim), all vertices [ps] | 62.1 | 62.2 | 62.4 |
| RMS, signal vertex [ps] | 20.0 | 20.3 | 20.5 |
| fraction within 30 ps | 0.473 | 0.471 | 0.471 |
| vertices with a time (1000 events) | 109962 | 110031 | 110250 |

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

Results, 300 events, 1 thread, Intel Xeon W-2155 (turbo on, governor powersave), input on a local HDD:

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
