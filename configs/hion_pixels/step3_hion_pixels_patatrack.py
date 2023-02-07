# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step3 -s RAW2DIGI,L1Reco,RECO,PAT,VALIDATION:@standardValidationNoHLT+@miniAODValidation,DQM:@standardDQMFakeHLT+@miniAODDQM --conditions auto:phase1_2021_realistic_hi -n 2 --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --eventcontent RECOSIM,MINIAODSIM,DQM --era Run3_pp_on_PbPb --procModifiers genJetSubEvent --filein file:step2.root --fileout file:step3.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_pp_on_PbPb_cff import Run3_pp_on_PbPb
from Configuration.ProcessModifiers.genJetSubEvent_cff import genJetSubEvent
from Configuration.ProcessModifiers.gpu_cff import gpu
from Configuration.ProcessModifiers.pixelNtupletFit_cff import pixelNtupletFit

process = cms.Process('RECO',Run3_pp_on_PbPb,genJetSubEvent,gpu,pixelNtupletFit)

from pbpb import pbpb
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')

options.register ('n',
                  100, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "n")
options.register ('threads',
                  8,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,
                  "threads")
options.register ('gpu',
                  True,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                  "gpu")
options.register ('timing',
                  False,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                  "timing")
options.parseArguments()


reconame = 'file:hionpixel'
dqmname = 'file:hionpixel_inDQM'

suff = "_" + str(options.n)  + "_t" + str(options.threads)

acc = "cpu"
if options.gpu:
	suff = suff + "_gpu"
	acc = "*"

reconame = reconame + suff + ".root"
dqmname  = dqmname  + suff + ".root"

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.n),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
 fileNames = cms.untracked.vstring(pbpb
),
 skipEvents = cms.untracked.uint32(0),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring(acc),
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
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:2'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(reconame),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('file:step3_inMINIAODSIM.root'),

    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVerticesWithBS__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(dqmname),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# process.hiPixelTracksSoA.cuda.trackQualityCuts.chi2MaxPt = cms.double(10)
# process.hiPixelTracksSoA.cuda.trackQualityCuts.chi2Coeff = cms.vdouble(0.9,1.8)
# process.hiPixelTracksSoA.cuda.trackQualityCuts.chi2Scale = cms.double(8)
# process.hiPixelTracksSoA.cuda.trackQualityCuts.tripletMinPt = cms.double(0)
# process.hiPixelTracksSoA.cuda.trackQualityCuts.tripletMaxTip = cms.double(0.1)
# process.hiPixelTracksSoA.cuda.trackQualityCuts.tripletMaxZip = cms.double(6)
# process.hiPixelTracksSoA.cuda.trackQualityCuts.quadrupletMinPt = cms.double(0)
# process.hiPixelTracksSoA.cuda.trackQualityCuts.quadrupletMaxTip = cms.double(0.5)
# process.hiPixelTracksSoA.cuda.trackQualityCuts.quadrupletMaxZip = cms.double(6)
#
# process.hiPixelTracksSoA.cpu.trackQualityCuts.chi2MaxPt = cms.double(10)
# process.hiPixelTracksSoA.cpu.trackQualityCuts.chi2Coeff = cms.vdouble(0.9,1.8)
# process.hiPixelTracksSoA.cpu.trackQualityCuts.chi2Scale = cms.double(8)
# process.hiPixelTracksSoA.cpu.trackQualityCuts.tripletMinPt = cms.double(0)
# process.hiPixelTracksSoA.cpu.trackQualityCuts.tripletMaxTip = cms.double(0.1)
# process.hiPixelTracksSoA.cpu.trackQualityCuts.tripletMaxZip = cms.double(6)
# process.hiPixelTracksSoA.cpu.trackQualityCuts.quadrupletMinPt = cms.double(0)
# process.hiPixelTracksSoA.cpu.trackQualityCuts.quadrupletMaxTip = cms.double(0.5)
# process.hiPixelTracksSoA.cpu.trackQualityCuts.quadrupletMaxZip = cms.double(6)


# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '125X_mcRun3_2022_realistic_HI_v3', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)

process.reconstruction_step = cms.Path(process.reconstruction)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

import validation

process.trackValidatorHIonPixelTrackingOnly = validation.trackValidatorHIonPixelTrackingOnly.clone()
process.trackingParticlePixelTrackAsssociationHIon = validation.trackingParticlePixelTrackAsssociationHIon.clone()
process.tpClusterProducerHIon = validation.tpClusterProducerHIon.clone(pixelClusterSrc = cms.InputTag("siPixelClustersPreSplitting"))
process.quickTrackAssociatorByHitsHIon = validation.quickTrackAssociatorByHitsHIon.clone()
process.VertexAssociatorByPositionAndTracksHIon = validation.VertexAssociatorByPositionAndTracksHIon.clone()

process.tracksValidation = cms.Sequence(process.trackValidatorHIonPixelTrackingOnly)
process.tracksValidationTruth = cms.Task(process.VertexAssociatorByPositionAndTracksHIon, process.quickTrackAssociatorByHitsHIon, process.tpClusterProducerHIon, process.trackingParticleNumberOfLayersProducer, process.trackingParticlePixelTrackAsssociationHIon)

process.validation = cms.Sequence(process.tracksValidation,process.tracksValidationTruth)
process.validation_step = cms.EndPath(process.validation)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.validation_step,process.DQMoutput_step)
if options.timing:
    suff = suff + "_timing"
    process.consumer = cms.EDAnalyzer("GenericConsumer", eventProducts = cms.untracked.vstring("hiConformalPixelTracks"))
    process.consume_step = cms.EndPath(process.consumer)
    process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.consume_step)

# customisation of the process.

process.load( "HLTrigger.Timer.FastTimerService_cfi" )
# print a text summary at the end of the job
process.FastTimerService.printEventSummary         = False
process.FastTimerService.printRunSummary           = False
process.FastTimerService.printJobSummary           = True

# enable DQM plots
process.FastTimerService.enableDQM                 = True

# enable per-path DQM plots (starting with CMSSW 9.2.3-patch2)
process.FastTimerService.enableDQMbyPath           = True

# enable per-module DQM plots
process.FastTimerService.enableDQMbyModule         = True

# enable per-event DQM plots vs lumisection
process.FastTimerService.enableDQMbyLumiSection    = True
process.FastTimerService.dqmLumiSectionsRange      = 2500

# set the time resolution of the DQM plots
tr = 10000000000.
tp = 10000000000.
tm = 2000000000.
process.FastTimerService.dqmTimeRange              = tr
process.FastTimerService.dqmTimeResolution         = tr/100.0
process.FastTimerService.dqmPathTimeRange          = tp
process.FastTimerService.dqmPathTimeResolution     = tp/100.0
process.FastTimerService.dqmModuleTimeRange        = tm
process.FastTimerService.dqmModuleTimeResolution   = tm/100.0

# set the base DQM folder for the plots
process.FastTimerService.dqmPath                   = 'HLT/TimerService'
process.FastTimerService.enableDQMbyProcesses      = True


process.FastTimerService.dqmMemoryRange            = 1000000
process.FastTimerService.dqmMemoryResolution       =    5000
process.FastTimerService.dqmPathMemoryRange        = 1000000
process.FastTimerService.dqmPathMemoryResolution   =    5000
process.FastTimerService.dqmModuleMemoryRange      =  100000
process.FastTimerService.dqmModuleMemoryResolution =     500

process.FastTimerService.writeJSONSummary = cms.untracked.bool(True)

wf_name = 'hionpixel_' + suff + '.json'

process.FastTimerService.jsonFileName = cms.untracked.string(wf_name)

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions

# customisation of the process.
#process.hiPixelTracksCUDA.doClusterCut = False
#process.hiPixelTracksSoA.cpu.doClusterCut = False
# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# End of customisation functions

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
