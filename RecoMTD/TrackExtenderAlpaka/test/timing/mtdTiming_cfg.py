# Timing of the full reconstruction with the legacy or the Alpaka MTD chain.
#
# Same reconstruction as ../validation/mtdValidation_cfg.py (same options: timing=, vertices=, threads=,
# maxEvents=, inputFiles=, output=), with:
#   - no validation, no DQM, no output file;
#   - a GenericConsumer that asks for the RECOSIM products made in this job, so that the same modules
#     run as in a RECO step (the input products are not asked for, so nothing extra is read from disk);
#   - FastTimerService (time per module, json summary <output>.json) and ThroughputService;
#   - threads = streams.
import os
import sys

import FWCore.ParameterSet.Config as cms

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'validation'))
from mtdValidation_cfg import process, options  # noqa: E402

# no validation, no DQM output
for name in ['validation_step', 'DQMoutput_step']:
    process.schedule.remove(getattr(process, name))
    delattr(process, name)
del process.DQMoutput

# ask for every RECOSIM product made by this process (keep statements only, process name RERECO)
keeps = [c.split()[1] for c in process.RECOSIMEventContent.outputCommands if c.startswith('keep')]
patterns = sorted({'_'.join(k.split('_')[:3] + [process.name_()]) for k in keeps if k.count('_') == 3})
process.recoConsumer = cms.EDAnalyzer('GenericConsumer', eventProducts=cms.untracked.vstring(patterns))
process.consumer_step = cms.EndPath(process.recoConsumer)
process.schedule.append(process.consumer_step)

process.options.numberOfStreams = options.threads
process.options.wantSummary = False
process.MessageLogger.cerr.FwkReport.reportEvery = 10
# the timing summaries are INFO messages, hidden by default
for category in ['FastReport', 'ThroughputService']:
    setattr(process.MessageLogger.cerr, category, cms.untracked.PSet(limit=cms.untracked.int32(10000000)))

process.load('HLTrigger.Timer.FastTimerService_cfi')
process.FastTimerService.enableDQM = False
process.FastTimerService.printRunSummary = False
process.FastTimerService.printJobSummary = True
process.FastTimerService.writeJSONSummary = True
process.FastTimerService.jsonFileName = options.output.replace('.root', '') + '.json'

process.ThroughputService = cms.Service('ThroughputService',
    eventRange=cms.untracked.uint32(100000),
    eventResolution=cms.untracked.uint32(1),
    eventSkip=cms.untracked.uint32(int(os.environ.get('MTD_TIMING_SKIP', '5'))),
    printEventSummary=cms.untracked.bool(True),
    enableDQM=cms.untracked.bool(False))
