import FWCore.ParameterSet.Config as cms

def customiseScoutngWithFullTracks(process):

    process.hltPixelTracksZetaCleanLowPt = cms.EDProducer( "TrackWithVertexSelector",
        normalizedChi2 = cms.double( 999999.0 ),
        numberOfValidHits = cms.uint32( 0 ),
        zetaVtx = cms.double( 0.3 ),
        etaMin = cms.double( 0.0 ),
        rhoVtx = cms.double( 0.2 ),
        ptErrorCut = cms.double( 5.0 ),
        dzMax = cms.double( 999.0 ),
        etaMax = cms.double( 5.0 ),
        quality = cms.string( "any" ),
        copyTrajectories = cms.untracked.bool( False ),
        nSigmaDtVertex = cms.double( 0.0 ),
        timesTag = cms.InputTag( "" ),
        ptMin = cms.double( -1.0 ),
        ptMax = cms.double( 6.0 ),
        d0Max = cms.double( 999.0 ),
        copyExtras = cms.untracked.bool( False ),
        nVertices = cms.uint32( 2 ),
        vertexTag = cms.InputTag( "hltPixelVertices" ),
        src = cms.InputTag( "hltPixelTracks" ),
        vtxFallback = cms.bool( True ),
        numberOfLostHits = cms.uint32( 999 ),
        numberOfValidPixelHits = cms.uint32( 3 ),
        timeResosTag = cms.InputTag( "" ),
        useVtx = cms.bool( True )
    )

    process.hltPixelTracksZetaCleanHighPt = cms.EDProducer( "TrackWithVertexSelector",
        normalizedChi2 = cms.double( 999999.0 ),
        numberOfValidHits = cms.uint32( 0 ),
        zetaVtx = cms.double( 0.3 ),
        etaMin = cms.double( 0.0 ),
        rhoVtx = cms.double( 0.2 ),
        ptErrorCut = cms.double( 5.0 ),
        dzMax = cms.double( 999.0 ),
        etaMax = cms.double( 5.0 ),
        quality = cms.string( "any" ),
        copyTrajectories = cms.untracked.bool( False ),
        nSigmaDtVertex = cms.double( 0.0 ),
        timesTag = cms.InputTag( "" ),
        ptMin = cms.double( 6.0 ),
        ptMax = cms.double( 500.0 ),
        d0Max = cms.double( 999.0 ),
        copyExtras = cms.untracked.bool( False ),
        nVertices = cms.uint32( 2 ),
        vertexTag = cms.InputTag( "hltPixelVertices" ),
        src = cms.InputTag( "hltPixelTracks" ),
        vtxFallback = cms.bool( True ),
        numberOfLostHits = cms.uint32( 999 ),
        numberOfValidPixelHits = cms.uint32( 3 ),
        timeResosTag = cms.InputTag( "" ),
        useVtx = cms.bool( True )
    )

    process.hltPixelOnlyPFSeedsFromPixelTracks = cms.EDProducer("SeedGeneratorFromProtoTracksEDProducer",
        InputCollection = cms.InputTag("hltPixelTracksZetaCleanHighPt"),
        InputVertexCollection = cms.InputTag(""),
        SeedCreatorPSet = cms.PSet(
            refToPSet_ = cms.string('HLTSeedFromProtoTracks')
        ),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        includeFourthHit = cms.bool(True),
        originHalfLength = cms.double(1000.0),
        originRadius = cms.double(1000.0),
        useEventsWithNoVertex = cms.bool(True),
        usePV = cms.bool(False),
        useProtoTrackKinematics = cms.bool(False)
    )

    process.hltPixelOnlyPFCkfFullTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
        MeasurementTrackerEvent = cms.InputTag("hltSiStripClusters"),
        NavigationSchool = cms.string('SimpleNavigationSchool'),
        RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
        SimpleMagneticField = cms.string('ParabolicMf'),
        TrajectoryBuilder = cms.string(''),
        TrajectoryBuilderPSet = cms.PSet(
            refToPSet_ = cms.string('HLTIter0GroupedCkfTrajectoryBuilderIT')
        ),
        TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
        TransientInitialStateEstimatorParameters = cms.PSet(
            numberMeasurementsForFit = cms.int32(4),
            propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
            propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
        ),
        cleanTrajectoryAfterInOut = cms.bool(False),
        doSeedingRegionRebuilding = cms.bool(False),
        maxNSeeds = cms.uint32(100000),
        maxSeedsBeforeCleaning = cms.uint32(1000),
        reverseTrajectories = cms.bool(False),
        src = cms.InputTag("hltPixelOnlyPFSeedsFromPixelTracks"),
        useHitsSplitting = cms.bool(False)
    )

    process.hltPixelOnlyPFFullTracks = cms.EDProducer("TrackProducer",
        AlgorithmName = cms.string('hltIter0'),
        Fitter = cms.string('hltESPFittingSmootherIT'),
        GeometricInnerState = cms.bool(True),
        MeasurementTracker = cms.string(''),
        MeasurementTrackerEvent = cms.InputTag("hltSiStripClusters"),
        NavigationSchool = cms.string(''),
        Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
        SimpleMagneticField = cms.string('ParabolicMf'),
        TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
        TrajectoryInEvent = cms.bool(False),
        alias = cms.untracked.string('ctfWithMaterialTracks'),
        beamSpot = cms.InputTag("hltOnlineBeamSpot"),
        clusterRemovalInfo = cms.InputTag(""),
        src = cms.InputTag("hltPixelOnlyPFCkfFullTrackCandidates"),
        useHitsSplitting = cms.bool(False),
        useSimpleMF = cms.bool(True)
    )

    process.hltPixelOnlyPFMuonMerging = cms.EDProducer( "TrackListMerger",
        ShareFrac = cms.double( 0.19 ),
        writeOnlyTrkQuals = cms.bool( False ),
        MinPT = cms.double( 0.05 ),
        allowFirstHitShare = cms.bool( True ),
        copyExtras = cms.untracked.bool( True ),
        Epsilon = cms.double( -0.001 ),
        selectedTrackQuals = cms.VInputTag( 'hltIterL3MuonTracks','hltPixelTracksZetaCleanLowPt','hltPixelOnlyPFFullTracks' ),
        indivShareFrac = cms.vdouble( 1.0, 1.0, 1.0 ),
        MaxNormalizedChisq = cms.double( 1000.0 ),
        copyMVA = cms.bool( False ),
        FoundHitBonus = cms.double( 5.0 ),
        LostHitPenalty = cms.double( 20.0 ),
        setsToMerge = cms.VPSet(
          cms.PSet(  pQual = cms.bool( False ),
            tLists = cms.vint32( 0, 1, 2 )
          )
        ),
        MinFound = cms.int32( 3 ),
        hasSelector = cms.vint32( 0, 0, 0 ),
        TrackProducers = cms.VInputTag( 'hltIterL3MuonTracks','hltPixelTracksZetaCleanLowPt','hltPixelOnlyPFFullTracks'),
        trackAlgoPriorityOrder = cms.string( "hltESPTrackAlgoPriorityOrder" ),
        newQuality = cms.string( "confirmed" )
    )


    process.HLTTrackReconstructionForPixelOnlyPF = cms.Sequence(
                          process.HLTDoLocalPixelSequence +
                          process.HLTRecopixelvertexingSequence +
                          process.hltPixelTracksZetaCleanLowPt + #new
                          process.hltPixelTracksZetaCleanHighPt + #new
                          process.hltPixelOnlyPFSeedsFromPixelTracks + #new
                          process.hltPixelOnlyPFCkfFullTrackCandidates + #new
                          process.hltPixelOnlyPFFullTracks + #new
                          process.hltPixelOnlyPFMuonMerging +
                          process.hltPixelOnlyMuonLinks +
                          process.hltPixelOnlyMuons )

    return process
