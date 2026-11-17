import FWCore.ParameterSet.Config as cms


def _replaceByAlias(process, label, alias):
    """an EDAlias cannot replace a module that is still in a task, sequence or path: take it out first"""
    module = getattr(process, label)
    for holder in list(process.tasks_().values()) + list(process.sequences_().values()) + \
            list(process.paths_().values()) + list(process.endpaths_().values()):
        while holder.remove(module):
            pass
    delattr(process, label)
    setattr(process, label, alias)


def customizeMtdTimingAlpaka(process, trackClones=True, deviceMVA=True):
    """Replace the legacy TrackExtenderWithMTD by the Alpaka MTD timing (no refit).

    mtdEventInput packs everything the device chain needs from the event into one host product, so
    that it goes to the device with a single copy; the BTL geometry is an event setup product (copied
    once per geometry). The new output module keeps the legacy label 'trackExtenderWithMTD' and the
    legacy product names, so every downstream module (mtdTrackQualityMVA, tofPID*, 4D vertices, TICL,
    MTD validation) reads the new products without any change.

    deviceMVA: the MTD track-quality BDT is computed on the device too, and 'mtdTrackQualityMVA'
    becomes an alias of its output (trackExtenderWithMTD:mtdQualMVA).
    """
    from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesAlpaka_cfi import pvTracksSoA as _pvTracksSoA
    process.mtdBtlGeometryESProducer = cms.ESProducer('MtdBtlGeometryESProducer@alpaka')
    process.mtdEventInput = cms.EDProducer('MtdEventInputProducer@alpaka',
        tracksSrc=cms.InputTag('generalTracks'),
        hitsSrc=cms.InputTag('mtdTrackingRecHits'),
        beamSpotSrc=cms.InputTag('offlineBeamSpot'),
        # tracks for the device vertexing: same selection as the standard producer
        TkFilterParameters=_pvTracksSoA.TkFilterParameters.clone())
    process.mtdTimingProducer = cms.EDProducer('MtdTimingProducer@alpaka',
        inputSrc=cms.InputTag('mtdEventInput'))
    newExtender = cms.EDProducer('MtdTimingToValueMaps',
        tracksSrc=cms.InputTag('generalTracks'),
        timingSrc=cms.InputTag('mtdTimingProducer'),
        trackClones=cms.bool(trackClones),
        hitsSrc=cms.InputTag('mtdTrackingRecHits'))
    process.globalReplace('trackExtenderWithMTD', newExtender)
    process.mtdTimingTask = cms.Task(process.mtdEventInput, process.mtdTimingProducer)
    if deviceMVA and hasattr(process, 'mtdTrackQualityMVA'):
        _replaceByAlias(process, 'mtdTrackQualityMVA', cms.EDAlias(trackExtenderWithMTD=cms.VPSet(
            cms.PSet(type=cms.string('floatedmValueMap'), fromProductInstance=cms.string('mtdQualMVA'),
                     toProductInstance=cms.string('mtdQualMVA')))))
    # the modules run when their products are consumed (unscheduled), like the rest of the reconstruction
    if process.schedule is not None:
        process.schedule.associate(process.mtdTimingTask)
    else:
        for path in process.paths.values():
            path.associate(process.mtdTimingTask)
    return process


def customizeMtdDeviceVertices(process):
    """3D vertices and the TOF PID on the device, on top of customizeMtdTimingAlpaka (call it first).

    unsortedOfflinePrimaryVertices: DA in blocks on the device, from the packed event input. One clustering,
    two fits (no beam spot "" and with beam spot "WithBS"), vertex time from the MTD times of the tracks
    (fromTracksPID) on the device. The two collections come back as reco::Vertex and keep the legacy label
    through an alias, so the sorters, the 4D producer and the rest read them unchanged.

    tofPID4DnoPID: TOFPIDProducer on the device, on the device 3D vertices; same parameters as the legacy
    module, same product names.
    """
    import RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesAlpaka_cfi as _pv
    process.pvVertexSoA = _pv.pvVertexSoA.clone(
        TrackLabel='mtdEventInput',
        inputType='mtdEventInput',
        alsoWithBS=cms.bool(True),
        timingSrc=cms.InputTag('mtdTimingProducer'))
    process.unsortedOfflinePrimaryVerticesFromSoA = _pv.unsortedOfflinePrimaryVerticesFromSoA.clone(
        soaVertex='pvVertexSoA', vertexTimeParameters=dict(algorithm='fromSoA'))
    process.unsortedOfflinePrimaryVerticesFromSoAWithBS = _pv.unsortedOfflinePrimaryVerticesFromSoA.clone(
        soaVertex=cms.InputTag('pvVertexSoA', 'WithBS'), vertexTimeParameters=dict(algorithm='fromSoA'))
    _replaceByAlias(process, 'unsortedOfflinePrimaryVertices', cms.EDAlias(
        unsortedOfflinePrimaryVerticesFromSoA=cms.VPSet(cms.PSet(
            type=cms.string('recoVertexs'), fromProductInstance=cms.string(''), toProductInstance=cms.string(''))),
        unsortedOfflinePrimaryVerticesFromSoAWithBS=cms.VPSet(cms.PSet(
            type=cms.string('recoVertexs'), fromProductInstance=cms.string(''),
            toProductInstance=cms.string('WithBS')))))

    legacy = process.tofPID4DnoPID
    process.tofPidProducer = cms.EDProducer('TofPidProducer@alpaka',
        inputSrc=cms.InputTag('mtdEventInput'),
        timingSrc=cms.InputTag('mtdTimingProducer'),
        vtxsSrc=cms.InputTag('pvVertexSoA'),
        **{p: getattr(legacy, p) for p in ['vtxMaxSigmaT', 'maxDz', 'maxDtSignificance', 'minProbHeavy',
                                           'fixedT0Error', 'probPion', 'probKaon', 'probProton',
                                           'vertexReassignment', 'MVASel', 'minTrackTimeQuality']})
    process.globalReplace('tofPID4DnoPID', cms.EDProducer('TofPidToValueMaps',
        tracksSrc=cms.InputTag('generalTracks'), tofPidSrc=cms.InputTag('tofPidProducer')))

    process.mtdDeviceVertexTask = cms.Task(process.pvVertexSoA, process.unsortedOfflinePrimaryVerticesFromSoA,
                                           process.unsortedOfflinePrimaryVerticesFromSoAWithBS,
                                           process.tofPidProducer)
    if process.schedule is not None:
        process.schedule.associate(process.mtdDeviceVertexTask)
    else:
        for path in process.paths.values():
            path.associate(process.mtdDeviceVertexTask)
    return process
