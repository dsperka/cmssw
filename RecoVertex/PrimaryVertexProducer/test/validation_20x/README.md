# DA in blocks (C++ and Alpaka) and the GNN vertex slots on CMSSW 20_1_X: validation

This branch puts three things on one base and compares them with the same validation:

1. **Base:** `CMSSW_20_1_X_2026-11-02-1100` (integration build, `el9_amd64_gcc14`).
2. **PR 51862** (vertex-slot GNN, behind the `vertexSlotGNN` process modifier), applied unchanged.
3. **DA in blocks** (PR 50766 and the developments of the 16_1 campaign): the C++ in-blocks clusterizer
   with a z-error-widened merge threshold, and the Alpaka CPU/GPU port (`PortableTrackSoAProducer`,
   `PrimaryVertexProducerPortable`, `RecoVertexFromSoAProducer`), with the review changes applied.

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
   with the GNN modifier on, plus the C++ in-blocks and Alpaka chains. Writes the FEVT file (about 170 MB per event: 34 GB for 200 events).
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
The first 100 events alone (`results/summary_per_event.txt`, one file) give the same picture to within statistical noise.

| | D.A. 3D | D.A. 4D | C++ In Blocks | Alpaka CPU | GNN |
|---|---|---|---|---|---|
| reco, all | 138.5 | 137.4 | 143.6 | 141.6 | 154.4 |
| pass ndof | 129.3 | 120.1 | 136.0 | 126.0 | 112.6 |
| real | 119.2 | 111.7 | 120.5 | 116.1 | 98.5 |
| fake | 10.1 | 8.4 | 15.5 | 9.9 | 14.2 |
| split | 8.4 | 6.2 | 12.2 | 6.6 | 12.7 |
| real / (real + fake) | 0.922 | 0.930 | 0.886 | 0.921 | 0.874 |
| split fraction | 0.065 | 0.052 | 0.090 | 0.052 | 0.113 |
| time, ms / event | 779 | 838 | 230 | 314 | 1343 |

Time is the sum of the producer modules of each chain from the reconstruction log (including the vertex sorter, about 110 ms in each),
single run, one thread.

What it shows:
- C++ in blocks is 3.4 times faster than the legacy 3D DA; Alpaka CPU on the serial backend is 2.5 times faster, with purity comparable to the 3D DA.
- **C++ in blocks and Alpaka CPU disagree here** (fake 15.5 vs 9.9 per event, split 12.2 vs 6.6). They agreed closely in the 16_1 campaign
  under the position-only matching of `PrimaryVertexAnalyzer4PUSlimmed`. A likely cause is the different fitters (weighted-mean in C++, the device fitter in Alpaka);
  this was not tested.
- The GNN is below the DA variants on real vertices and purity on this sample (D128, PU150). 42 of its 154 vertices per event fail the ndof cut (200 events: 154.4 reco, 112.6 pass).
  The PR targets D121 PU200, and the model may be trained for that.

## Files
- `results/summary_per_event_200events.txt`, `results/DQM_6way_harvested_200events.root`: the numbers above and the harvested file (200 events).
  The `*_200events` and plain versions of each file are the 200-event and first-100-event results.
- `galleries/pv_dqm_gallery_6way_200events.html`: all plots shared by the five collections (self-contained HTML; the plain `pv_dqm_gallery_6way.html` is the 100-event one). `phone/`: phone-sized pages.
- `render_dqm_gallery.py`, `dqm_labels.py`: the gallery generator. Options added for this work: `--base-dir`, `--display KEY=NAME`, a fifth series (`--extra`).
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
Alpaka GPU (needs the lxplus recipe above); the C++ versus Alpaka difference; a D121 PU200 sample for the GNN; `PrimaryVertexAnalyzer4PUSlimmed` rebinning (vertex-count histograms were empty with the default range for PU150; not needed for the MTD validation).
