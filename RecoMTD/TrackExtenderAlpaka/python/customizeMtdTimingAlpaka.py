import FWCore.ParameterSet.Config as cms


def customizeMtdTimingAlpaka(process, trackClones=True):
    """Replace the legacy TrackExtenderWithMTD by the Alpaka MTD timing (no refit).

    The new output module keeps the legacy label 'trackExtenderWithMTD' and the legacy product
    names, so every downstream module (mtdTrackQualityMVA, tofPID*, 4D vertices, TICL, MTD
    validation) reads the new products without any change.
    """
    process.mtdTimingProducer = cms.EDProducer('MtdTimingProducer@alpaka',
        tracksSrc=cms.InputTag('generalTracks'),
        hitsSrc=cms.InputTag('mtdTrackingRecHits'))
    newExtender = cms.EDProducer('MtdTimingToValueMaps',
        tracksSrc=cms.InputTag('generalTracks'),
        timingSrc=cms.InputTag('mtdTimingProducer'),
        trackClones=cms.bool(trackClones),
        hitsSrc=cms.InputTag('mtdTrackingRecHits'))
    process.globalReplace('trackExtenderWithMTD', newExtender)
    process.mtdTimingTask = cms.Task(process.mtdTimingProducer)
    # the module runs when its product is consumed (unscheduled), like the rest of the reconstruction
    if process.schedule is not None:
        process.schedule.associate(process.mtdTimingTask)
    else:
        for path in process.paths.values():
            path.associate(process.mtdTimingTask)
    return process
