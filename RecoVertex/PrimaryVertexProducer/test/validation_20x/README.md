# DA in blocks (C++ and Alpaka) and the GNN vertex slots on CMSSW 20_1_X: validation

This branch puts three things on one base and compares them with the same validation:

1. **Base:** `CMSSW_20_1_X_2026-11-02-1100` (integration build, `el9_amd64_gcc14`).
2. **PR 51862** (vertex-slot GNN, behind the `vertexSlotGNN` process modifier), applied unchanged.
3. **DA in blocks** (PR 50766 and the developments of the 16_1 campaign): the C++ in-blocks clusterizer
   with a z-error-widened merge threshold, and the Alpaka CPU/GPU port (`PortableTrackSoAProducer`,
   `PrimaryVertexProducerPortable`, `RecoVertexFromSoAProducer`, which also computes the vertex time), with the review changes applied.

Commits on the branch, in order:
- `PR 51862 as submitted` (the 40 files of the PR, no change).
- `Port DA in blocks`: our developments on top, with the conflicts resolved (see below).
- `validation_20x`: the macros, instructions and results in this directory.

## What had to be merged by hand
Both PRs edit the same files in `DataFormats/VertexSoA` and `RecoVertex/PrimaryVertexProducer/plugins/BuildFile.xml`.
- Dictionaries (`classes.h`, `classes_def.xml`, `src/alpaka/classes_{cuda,rocm}{.h,_def.xml}`): both sets of entries kept.
- Serialisation: the Alpaka plugin uses the new `DEFINE_TRIVIAL_SERIALISER_PORTABLE_PLUGIN(host, device)` macro,
  for the vertexing collections as well as the GNN ones.
- `plugins/BuildFile.xml`: the structure of PR 51862 is kept; our dependencies and
  `<flags CUDA_FLAGS="-Xptxas -maxrregcount=100"/>` (register-pressure control for the device kernels) are added.
- `mintrkweight` is a parameter of the clusterizer, with default 0.5 (the previously hardcoded value). Setting another value from
  python needs a typed value: `mintrkweight = cms.double(0.25)`.
The tree built with no errors (serial, CUDA and ROCm plugins). No Alpaka API changes were needed between 16_1 and this base.

## The model
The PR needs `RecoVertex/PrimaryVertexProducer/data/vertexSlotGNN.pt`, which comes from cms-data and is **not** committed here.
Copy it to that path before running. It was checked to accept 13 features and return 5 tensors with 220 slots.

## Environment
```
export SCRAM_ARCH=el9_amd64_gcc14
source /cvmfs/cms-ib.cern.ch/sw/x86_64/week1/cmsset_default.sh
scram p CMSSW CMSSW_20_1_X_2026-11-02-1100 && cd CMSSW_20_1_X_2026-11-02-1100/src && cmsenv
git cms-addpkg ...      # or check this branch out over src/ (only the packages touched by the branch are needed)
scram b -j 16
```
**OpenSSL.** Every `cmsRun` of this base aborts at start-up on a machine with OpenSSL older than 3.4 (AlmaLinux 9.3 has 3.0.7):
the xrootd plugin needs `OPENSSL_3.4.0`. Workaround used here, with no change to the system: install an OpenSSL 3.x >= 3.4
(for example with micromamba) into your own area, symlink only `libssl.so.3` and `libcrypto.so.3` into a directory, and put it first in `LD_LIBRARY_PATH`.

## Run
Input used: `RelValTTbar_14TeV/CMSSW_20_1_0_pre3-PU_150X_mcRun4_realistic_v1_STD_D128_RegeneratedGS_PU_20260903_192421_RV19-v1/GEN-SIM-RECO`,
files `01313427-f364-4b67-9190-ea6bafe583b6.root` and `044bab19-4f2d-4f18-842b-32d35b81e270.root` (17 GB and 100 events each), 200 events in all.
They were copied locally, as there was no network access to the grid.
```
./run_validation.sh /path/to/01313427-f364-4b67-9190-ea6bafe583b6.root,/path/to/044bab19-4f2d-4f18-842b-32d35b81e270.root 1000 cpu
python3 summarize_dqm.py DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root
./make_gallery_6way.sh DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root gallery_6way
```
`run_validation.sh` runs three steps in the current directory:
1. `reco_cfg_template.py`: the PR's reconstruction config adapted to D128 (era `Phase2C26I13M9`, `GeometryExtendedRun4D128Reco`),
   with the GNN modifier on, plus the C++ in-blocks and Alpaka chains. Writes the FEVT file with a slim event content (only what the validation reads, about 66 MB per event; the full `FEVTDEBUGHLT` content is 170 MB per event; the slim file gives the same validation histograms and the job runs about 1.7 times faster with 6 threads, see the comment in the template).
2. `mtdValidation_6way_cfg.py`: the PR's `mtdValidation_3way_cfg.py` extended by two modules for the new collections.
3. `mtdHarvesting_6way_cfg.py`.
The first argument can be a local path, a `/store/...` LFN or a `root://` URL. Third argument: `cpu` (default) or `gpu-nvidia`.
The GPU was not run for the results below: this base uses CUDA 13.x, which does not support the pre-Turing GPUs (TITAN V, Quadro P400) of the machine used. It needs a Turing-or-newer GPU, see the lxplus recipe below.

Collections produced and validated:
| folder in the output | collection | how |
|---|---|---|
| `3D` | `offlinePrimaryVertices` | standard DA |
| `4D` | `offlinePrimaryVertices4D` | standard 4D DA |
| `GNN` | `offlinePrimaryVerticesGNN` | PR 51862 |
| `InBlocks` | `offlinePrimaryVerticesInBlocks` | standard producer with `runInBlocks`, 512-track blocks, overlap 0.5, `mintrkweight` 0.25, weighted-mean fitter |
| `Alpaka` | `offlinePrimaryVerticesAlpaka` | `tracksSoA` -> `vertexSoA` -> `RecoVertexFromSoAProducer`, then the standard sorter |

## Results on this sample (200 events, CPU, serial Alpaka backend)
Per-event averages from `Primary4DVertexValidation` (track-based matching, ndof cut). 189.4 simulated vertices per event.
Alpaka runs with `maxEta = 4.0`, like the C++ in blocks, so both cover the same tracks.

| | D.A. 3D | D.A. 4D | C++ In Blocks | Alpaka CPU | GNN |
|---|---|---|---|---|---|
| reco, all | 138.5 | 137.4 | 143.4 | 143.4 | 154.4 |
| pass ndof | 129.3 | 120.1 | 135.9 | 135.9 | 112.6 |
| real | 119.2 | 111.7 | 120.5 | 120.5 | 98.5 |
| fake | 10.1 | 8.4 | 15.4 | 15.4 | 14.2 |
| split | 8.4 | 6.2 | 12.1 | 12.0 | 12.7 |
| real / (real + fake) | 0.922 | 0.930 | 0.887 | 0.887 | 0.874 |
| split fraction | 0.065 | 0.052 | 0.089 | 0.088 | 0.113 |
| leading vertex, pileup track fraction | 0.288 | 0.250 | 0.302 | 0.301 | 0.208 |
| vertex time: RMS of t(rec) - t(sim), ns | 0.0568 | 0.0510 | 0.0587 | 0.0587 | 0.0516 |
| vertex time: fraction within 30 ps | 0.693 | 0.723 | 0.680 | 0.681 | 0.719 |
| time, ms / event | 752 | 834 | 237 | 366 | 1354 |

Time is the sum of the producer modules of each chain from the reconstruction log (including the vertex sorter, about 110 ms in each),
single run, one thread, no other jobs on the machine. The in-blocks chains include the vertex time fit.

What it shows:
- C++ in blocks is 3.3 times faster than the legacy 3D DA; Alpaka CPU on the serial backend is about 2 times faster.
- **C++ in blocks and Alpaka CPU agree.** Vertex counts per event are the same, the leading-vertex fractions differ by 0.0005, and the vertex track lists overlap by 99.5% on average
  (vertices per event 135.87 and 135.90 passing ndof, in the direct collection comparison). Earlier versions of the Alpaka code differed (fake 9.9 against 15.5 per event); the causes were
  (1) the arbitration resort read stale entries of the per-block order lists (duplicate vertices), (2) the vertex list compaction in `finalizeVertices` did not update `nV`
  (stale copies of valid vertices were fitted again), (3) the cross-block merge did not chain neighbouring candidates like the C++ code,
  (4) in thermalize, an extra update at convergence and a different reassignment criterion, (5) the split acceptance check, (6) float versus double storage of the track and vertex quantities,
  and (7) the C++ outlier rejection did not store the track weights (rejected tracks kept weight 1). The pileup fraction of the leading vertex is the one to compare: 0.302 against 0.301.
- **Vertex time.** All collections get the vertex time from the same algorithm (`vertexTimeParameters.algorithm = fromTracksPID`: an iterative fit over the pion, kaon and proton
  hypotheses, weighted by the vertex track weights). For C++ in blocks it is set in the `vertexCollections` of the clone (replacing the VPSet drops the setting that the `phase2_timing_layer`
  modifier adds to the standard producer); for Alpaka, `RecoVertexFromSoAProducer` takes the same `vertexTimeParameters` and runs the algorithm on each converted vertex. An earlier version of
  this validation had no vertex time in the in-blocks collections (RMS 0.107 ns, 19% within 30 ps). The in-blocks time resolution (0.0587 ns) is close to the 3D one (0.0568 ns); the 4D and GNN collections are better (0.051 ns).
- Both in-blocks versions have more fake and split vertices than the 3D and 4D DA, and fewer tracks left in no vertex (0.055 against 0.068 and 0.100).
- The GNN is below the DA variants on real vertices and purity on this sample (D128, PU150). 42 of its 154 vertices per event fail the ndof cut (200 events: 154.4 reco, 112.6 pass).
  The PR targets D121 PU200, and the model may be trained for that.

## Results on 1000 events (ten files, CPU, serial Alpaka backend)
Input: the first ten files of the same dataset (1000 events), Alpaka at its defaults (`mintrkweight` 0.25, `zmerge` 0.01, `uniquetrkweight` 0.8).
The sixth series of the gallery is Alpaka with `mintrkweight` 0.4 and `zmerge` 0.02 (a tuning tried in the parameter scan; it trades 3 leading vertices for fewer fake and split vertices).

| | D.A. 3D | D.A. 4D | GNN | C++ In Blocks | Alpaka CPU | Alpaka, mtw 0.4 zm 0.02 |
|---|---|---|---|---|---|---|
| real | 118.91 | 111.44 | 97.91 | 120.18 | 120.17 | 117.07 |
| fake | 9.83 | 8.33 | 14.11 | 15.32 | 15.26 | 13.57 |
| split | 8.20 | 6.24 | 12.65 | 11.98 | 11.91 | 10.17 |
| real / (real + fake) | 0.924 | 0.930 | 0.874 | 0.887 | 0.887 | 0.896 |
| leading vertex found (of 1000) | 959 | 938 | 884 | 949 | 949 | 946 |

C++ in blocks and Alpaka agree again (the same 949 leading vertices); the legacy 3D DA finds 10 more, partly a statistical fluctuation.
Neither the clusterizer nor the fitter alone accounts for the difference to the legacy DA (swap test on 200 events).
These numbers and the galleries were produced with the compact vertex-track-list data format (`VertexTrackSoA`; Alpaka CPU differs from the previous run by less than 0.05 vertices per event, as expected from the layout sensitivity of the CPU -O3 build).
Files: `results/summary_per_event_1000events.txt`, `results/DQM_6way_harvested_1000events.root`, `galleries/pv_dqm_gallery_6way_1000events.html`, `phone/six_collections_1000events.html`.
The same job with 6 threads gives histograms identical to the single-thread job (differences < 1e-15) and a reconstruction 1.7 times faster.

## Files
- `results/summary_per_event_200events.txt`, `results/DQM_6way_harvested_200events.root`: the numbers above and the harvested file (200 events).
- `galleries/pv_dqm_gallery_6way_200events.html`: all plots shared by the five collections (self-contained HTML). `phone/`: phone-sized pages.
- `render_dqm_gallery.py`, `dqm_labels.py`, `mtd_labels.py`: the gallery generator. `mtd_labels.py` gives every MTD validation histogram a plain-language title, an explanatory x-axis and y-axis title
  and a short description (what is plotted, which vertices and tracks, what the matching does); the original histogram name stays on each card as the subtitle.
  Options added for this work: `--base-dir`, `--display KEY=NAME`, up to seven series (`--extra`, `--extra2`, `--extra3`).
- `summarize_dqm.py`: the per-event table.

## Running on lxplus, with a GPU (gives the seven-way comparison)
```
ssh lxplus-gpu.cern.ch                                   # or an HTCondor job with request_gpus = 1
export SCRAM_ARCH=el9_amd64_gcc14
source /cvmfs/cms-ib.cern.ch/sw/x86_64/week1/cmsset_default.sh
scram list CMSSW_20_1_X_2026-11-02                       # the area is under one of the nweek-* directories; `scram p` finds it
scram p CMSSW CMSSW_20_1_X_2026-11-02-1100
cd CMSSW_20_1_X_2026-11-02-1100/src && cmsenv
git cms-init
git cms-merge-topic dsperka:from-CMSSW_20_1_X_2026-11-02-1100-da-in-blocks-and-gnn
mkdir -p RecoVertex/PrimaryVertexProducer/data && cp /path/to/vertexSlotGNN.pt RecoVertex/PrimaryVertexProducer/data/
scram b -j 16
cudaIsEnabled; echo "cudaIsEnabled rc=$?"             # must be 0 on the GPU node
```
Then, in two separate work directories with about 40 GB free each (for example on /eos):
```
V=$CMSSW_BASE/src/RecoVertex/PrimaryVertexProducer/test/validation_20x
D=/store/relval/CMSSW_20_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_D128_RegeneratedGS_PU_20260903_192421_RV19-v1/2830000
IN=$D/01313427-f364-4b67-9190-ea6bafe583b6.root,$D/044bab19-4f2d-4f18-842b-32d35b81e270.root
mkdir cpu gpu
(cd cpu && $V/run_validation.sh $IN 1000 cpu)
(cd gpu && $V/run_validation.sh $IN 1000 gpu-nvidia)
$V/make_gallery_7way.sh cpu/DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root gpu/DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root gallery_7way
python3 $V/summarize_dqm.py gpu/DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root
```
The GPU job runs everything again (the CPU-only collections too) with the Alpaka chain and the GNN on the device; the seven-way gallery uses its
Alpaka and GNN folders next to the full CPU job. Seven series: 3D DA, 4D DA, C++ in blocks, Alpaka CPU, Alpaka GPU, GNN CPU, GNN GPU.
No OpenSSL workaround should be needed there.

## Not done
Alpaka GPU on this release (needs the lxplus recipe above); a D121 PU200 sample for the GNN; `PrimaryVertexAnalyzer4PUSlimmed` rebinning (vertex-count histograms were empty with the default range for PU150; not needed for the MTD validation).
