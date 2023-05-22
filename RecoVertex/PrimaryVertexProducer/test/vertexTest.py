import FWCore.ParameterSet.Config as cms
from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import offlinePrimaryVertices, offlinePrimaryVerticesCUDA, offlinePrimaryVerticesCUDA_CPUFitter,  offlinePrimaryVerticesDumbFitter
#from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesCUDA_cfi import offlinePrimaryVertices as offlinePrimaryVerticesCUDA
#from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesCUDA_cfi import offlinePrimaryVertices as offlinePrimaryVerticesDumpFitter
import FWCore.ParameterSet.VarParsing as VarParsing

from HeterogeneousCore.CUDACore.SwitchProducerCUDA import SwitchProducerCUDA

process = cms.Process("Vertexing")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("HeterogeneousCore.CUDAServices.CUDAService_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('commons_cff')


options = VarParsing.VarParsing('analysis')

options.register ('n',
                  -1, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "n")
options.register ('threads',
                  1,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,
                  "threads")
options.register ('gpu',
                  False,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                  "gpu")
options.register ('timing',
                  False,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                  "timing")
options.register ('both',
                  False,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                  "gpuVScpu")
options.register ('verbose',
                  False,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                  "verbose")
options.parseArguments()

from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mcRun4_realistic_v3', '')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

suff = 'gpu'

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(options.threads),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)


if not options.gpu:
    process.options.accelerators = cms.untracked.vstring('cpu')
    suff = 'cpu'

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.n))

process.source = cms.Source("PoolSource",
fileNames = cms.untracked.vstring(
#'file:/data/user/cericeci/QCD_5/reco_1002.root',
#'/store/relval/CMSSW_12_4_0_pre3/RelValTTbar_14TeV/GEN-SIM-RECO/PU_123X_mcRun4_realistic_v11_2026D88PU200-v1/2580000/7781d089-b51a-495a-b1ba-384c15e90749.root',
#'/store/relval/CMSSW_12_4_0_pre3/RelValZMM_14_HI_2021/GEN-SIM-RECO/123X_mcRun3_2021_realistic_HI_v14-v1/2580000/633b39de-a6cf-45ca-9182-96be30f293fe.root',
#'/store/relval/CMSSW_12_4_0_pre3/RelValTTbar_14TeV/GEN-SIM-RECO/PU_123X_mcRun4_realistic_v11_2026D88PU200-v1/2580000/f6b68ca4-5b0e-42bb-b1d0-f94480067693.root',
#'/store/relval/CMSSW_12_4_0_pre3/RelValTTbar_14TeV/GEN-SIM-RECO/PU_123X_mcRun4_realistic_v11_2026D88PU200-v1/2580000/876a46e3-477e-4c53-8a4a-c16e7c8dee0b.root'
#'file:aca7b050-5990-4576-a9ee-f41ac82e5b86.root'
#'file:/eos/user/j/jshteren/PV/7781d089-b51a-495a-b1ba-384c15e90749.root'
#'file:/eos/user/c/cericeci/PV/QCD_5_update/reco_1002.root'
'file:/cms/data/store/user/dsperka/PV/RelValTTbar_14TeV/CMSSW_12_6_0_pre5/0a765d47-7b73-4fa2-bfc4-299186b6d952.root'
),
skipEvents=cms.untracked.uint32(0),
inputCommands = cms.untracked.vstring(
        'keep *','drop *_offlinePrimaryVertices_*_*','drop l1tTkPrimaryVertexs_L1TkPrimaryVertex__HLT'
    )
)

if options.timing:

    process.ThroughputService = cms.Service('ThroughputService',
        eventRange = cms.untracked.uint32(10),
        eventResolution = cms.untracked.uint32(1),
        printEventSummary = cms.untracked.bool(True),
        enableDQM = cms.untracked.bool(True),
        dqmPathByProcesses = cms.untracked.bool(False),
        dqmPath = cms.untracked.string('Throughput'),
        timeRange = cms.untracked.double(1000),
        timeResolution = cms.untracked.double(1)
    )

    process.MessageLogger.cerr.ThroughputService = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    )

    process.FastTimerService.writeJSONSummary = cms.untracked.bool(True)
    process.FastTimerService.jsonFileName = cms.untracked.string('resources_'+suff+'.json')

    if not options.gpu:
        process.IgProfService = cms.Service("IgProfService",
          reportFirstEvent            = cms.untracked.int32(0),
          reportEventInterval         = cms.untracked.int32(25),
          reportToFileAtPostEvent     = cms.untracked.string("| gzip -c > igdqm.%I.gz")
        )



if options.gpu:
    #process.vertex = offlinePrimaryVerticesDumbFitter.clone()
    process.vertex = offlinePrimaryVerticesCUDA.clone()
else:
    process.vertex = offlinePrimaryVerticesCUDA_CPUFitter.clone()

if options.verbose:
    process.vertex.verbose = cms.untracked.bool(True)

##############################CHANGED
process.vertex.vertexCollections = cms.VPSet(
     [cms.PSet(label=cms.string(""),
               algorithm=cms.string("WeightedMeanFitter"),
               chi2cutoff = cms.double(2.5),
               minNdof=cms.double(0.0),
               useBeamConstraint = cms.bool(False),
#               useBeamConstraint = cms.bool(True), # to be like the GPU version
               maxDistanceToBeam = cms.double(1.0)
               ),
      cms.PSet(label=cms.string("WithBS"),
               algorithm = cms.string('WeightedMeanFitter'),
               chi2cutoff = cms.double(2.5),
               minNdof=cms.double(2.0),
               useBeamConstraint = cms.bool(True),
               maxDistanceToBeam = cms.double(1.0),
               )
      ]
)
##############################CHANGED

process.vertexAnalysis.vertexRecoCollections  = cms.VInputTag("vertex")
process.pvMonitor.vertexLabel = cms.InputTag("vertex")

if options.both:
    suff = "gpuVScpu"

process.output = cms.OutputModule("PoolOutputModule",
    fileName= cms.untracked.string("file:test_"+suff+".root"),
    outputCommands = cms.untracked.vstring(
                                'drop *_*_*_*',
                                'keep *_demo_*_*'
    )
)

##DQM Output step
process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.output.fileName = 'file:test_'+suff+'.root'
process.DQMoutput.fileName = 'file:test_dqm_'+suff+'.root'

process.tracksValidationTruth = cms.Task(process.VertexAssociatorByPositionAndTracks, process.quickTrackAssociatorByHits, process.tpClusterProducer)
process.pvValidation = cms.Sequence(process.vertexAnalysis,process.tracksValidationTruth)
process.prevalidation_step = cms.Path(process.pvValidation)

process.DQMOfflineVertex = cms.Sequence(process.pvMonitor)
process.dqmoffline_step = cms.EndPath(process.DQMOfflineVertex)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

process.vertexing_step = cms.Path(process.vertex)
process.output_step = cms.EndPath(process.output)

process.schedule = cms.Schedule(process.vertexing_step,process.prevalidation_step,process.dqmoffline_step,process.DQMoutput_step,process.output_step)

if options.timing:

    process.consumer = cms.EDAnalyzer("GenericConsumer", eventProducts = cms.untracked.vstring("vertex"))
    process.consume_step = cms.EndPath(process.consumer)
    process.schedule.append(process.consume_step)

#else:
#    process.schedule = cms.Schedule(process.vertexing_step)

if options.both:

    process.vertex = offlinePrimaryVertices.clone()
    process.vertexCUDA = offlinePrimaryVerticesCUDA.clone()
    process.vertexing_step = cms.Path(process.vertex,process.vertexCUDA)
    process.consumerCPU = cms.EDAnalyzer("GenericConsumer", eventProducts = cms.untracked.vstring("vertex"))
    process.consumerGPU = cms.EDAnalyzer("GenericConsumer", eventProducts = cms.untracked.vstring("vertexCUDA"))
    process.consume_step = cms.EndPath(process.consumerCPU,process.consumerGPU)

    if not options.timing:
        process.vertexAnalysis.vertexRecoCollections = cms.VInputTag("vertex","vertexCUDA")

    process.schedule.append(process.consume_step)
process.vertexAnalysis.nPUbins = cms.uint32(130)
