import FWCore.ParameterSet.Config as cms

def addPixelLocalRecoSoAonCPU(process) :

    from RecoLocalTracker.SiPixelRecHits.siPixelRecHitSoAFromLegacy_cfi import siPixelRecHitSoAFromLegacy as siPixelRecHitsPreSplittingSoA

    process.siPixelRecHitsPreSplitting = siPixelRecHitsPreSplittingSoA.clone()
    process.siPixelRecHitsPreSplitting.isUpgrade = True
    process.siPixelRecHitsPreSplitting.convertToLegacy = True

    process.PixelCPEFastESProducer.isUpgrade = True

    if not hasattr(process,"siPixelClustersPreSplitting"):

        from RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizerPreSplitting_cfi import siPixelClustersPreSplitting
        process.siPixelClustersPreSplitting = siPixelClustersPreSplitting.clone()

    if not hasattr(process,"siPixelRecHitsPreSplitting") or not hasattr(process,"siPixelClustersPreSplitting"):

        process.consumer = cms.EDAnalyzer("GenericConsumer",eventProducts = cms.untracked.vstring('siPixelRecHitsPreSplitting'))

        if hasattr(process,'reconstruction_step'):
            process.reconstruction_step += process.siPixelClustersPreSplitting + process.siPixelRecHitsPreSplitting
        else:
                process.consume_step = cms.EndPath(process.consumer)
                process.schedule.append(process.consume_step)

    return process

def addPixelLocalRecoSoAonGPU(process) :

    process.load("HeterogeneousCore.CUDAServices.CUDAService_cfi")
    
    from RecoLocalTracker.SiPixelClusterizer.siPixelRawToClusterCUDA_cfi import siPixelRawToClusterCUDA
    from RecoLocalTracker.SiPixelRecHits.siPixelRecHitCUDA_cfi import siPixelRecHitCUDA as siPixelRecHitCUDA
    from EventFilter.SiPixelRawToDigi.siPixelDigisSoAFromCUDA_cfi import siPixelDigisSoAFromCUDA as siPixelDigisSoAFromCUDA
    from RecoLocalTracker.SiPixelClusterizer.siPixelDigisClustersFromSoA_cfi import siPixelDigisClustersFromSoA
    from RecoLocalTracker.SiPixelRecHits.siPixelRecHitFromCUDA_cfi import siPixelRecHitFromCUDA as siPixelRecHitFromCUDA

    process.siPixelClustersPreSplittingCUDA = siPixelRawToClusterCUDA.clone()
    process.siPixelClustersPreSplittingCUDA.isUpgrade = True
    process.siPixelClustersPreSplittingCUDA.isRun2 = False

    process.siPixelDigisSoA = siPixelDigisSoAFromCUDA.clone(
        src = "siPixelClustersPreSplittingCUDA",
    )

    process.siPixelClustersPreSplitting = siPixelDigisClustersFromSoA.clone()
    process.siPixelClustersPreSplitting.src = "siPixelDigisSoA"
    process.siPixelClustersPreSplitting.isUpgrade = True

    process.PixelCPEFastESProducer.isUpgrade = True

    process.siPixelRecHitsPreSplittingCUDA = siPixelRecHitCUDA.clone(
        beamSpot = "offlineBeamSpotToCUDA"
    )

    process.siPixelRecHitsPreSplitting = siPixelRecHitFromCUDA.clone();

    process.consumer = cms.EDAnalyzer("GenericConsumer",eventProducts = cms.untracked.vstring('siPixelRecHitsPreSplitting'))

    if hasattr(process,'reconstruction_step'):
        process.reconstruction_step += process.offlineBeamSpotToCUDA + process.siPixelClustersPreSplittingCUDA + process.siPixelDigisSoA + process.siPixelRecHitsPreSplittingCUDA
    else:
            process.consume_step = cms.EndPath(process.consumer)
            process.schedule.append(process.consume_step)


    return process
