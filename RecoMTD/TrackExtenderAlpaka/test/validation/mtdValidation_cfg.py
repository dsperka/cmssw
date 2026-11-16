# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 -s RAW2DIGI,L1Reco,RECO,RECOSIM,VALIDATION:mtdRecoValid --conditions auto:phase2_realistic_T35 --datatier DQMIO --eventcontent DQM --geometry ExtendedRun4D128 --era Phase2C26I13M9 --filein <relval GEN-SIM-RECO> -n 100 --nThreads 8 --no_exec
#
# MTD timing validation, legacy TrackExtenderWithMTD vs the Alpaka MTD timing module.
# Re-runs the full reconstruction on RelVal GEN-SIM-RECO (D128, PU200) and fills MtdTracksValidation and
# Primary4DVertexValidation. Options (cmsRun mtdValidation_cfg.py key=value):
#   timing=legacy|alpaka   which module fills trackExtenderWithMTD (default alpaka)
#   inputFiles=...         GEN-SIM-RECO files (comma separated), or inputFilesList=files.txt
#   maxEvents=N            default -1
#   output=out.root        DQMIO output
#   threads=N              default 8
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C26I13M9_cff import Phase2C26I13M9

process = cms.Process('RERECO',Phase2C26I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D128Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
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
    forceEventSetupCacheClearOnNewRun = cms.obsolete.untracked.bool,
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:mtdValidation_DQMIO.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.validation_step = cms.EndPath(process.mtdAssociationProducers + process.mtdTracksValid + process.vertices4DValid)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,process.validation_step,process.DQMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads = 8
process.options.numberOfStreams = 0

# customisation of the process.

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

# ---- MTD timing validation ----
# drop the products of the original RECO so that everything is re-reconstructed here
process.source.inputCommands = cms.untracked.vstring('keep *', 'drop *_*_*_RECO')
process.options.accelerators = ["cpu"]
process.MessageLogger.cerr.FwkReport.reportEvery = 50

# truth association needed by MtdTracksValidation (tracks <-> TrackingParticles by shared hits)
process.mtdValidTruthTask = cms.Task(process.tpClusterProducer, process.quickTrackAssociatorByHits,
                                     process.trackingParticleRecoTrackAsssociation)
process.validation_step.associate(process.mtdValidTruthTask)

# ---- command line options ----
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register('timing', 'alpaka', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 'legacy: TrackExtenderWithMTD with refit; alpaka: Alpaka MTD timing module (no refit)')
options.register('threads', 8, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'number of threads')
options.register('output', 'mtdValidation_DQMIO.root', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 'DQMIO output file (outputFile would get a _numEventN suffix)')
options.parseArguments()

process.source.fileNames = options.inputFiles
process.maxEvents.input = options.maxEvents
process.DQMoutput.fileName = 'file:' + options.output
process.options.numberOfThreads = options.threads

if options.timing == 'alpaka':
    from RecoMTD.TrackExtenderAlpaka.customizeMtdTimingAlpaka import customizeMtdTimingAlpaka
    process = customizeMtdTimingAlpaka(process)
elif options.timing != 'legacy':
    raise ValueError('timing must be legacy or alpaka, not %s' % options.timing)
