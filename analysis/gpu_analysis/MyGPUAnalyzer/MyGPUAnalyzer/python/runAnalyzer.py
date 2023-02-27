ouput_filename = 'dimuonSoA.root'

import FWCore.ParameterSet.Config as cms
process = cms.Process('DiMuonSoA')

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('Configuration.StandardSequences.Accelerators_cff') 
process.load('HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v4', '')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))

from data import files
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(files))

process.TFileService = cms.Service("TFileService",fileName = cms.string(ouput_filename))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(4)

process.selectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('slimmedMuons'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
                    ' && (pt > 2.)'
   ),
   filter = cms.bool(True)
)

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring('HLT_DoubleMu4_3_LowMass_v*'),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.producer = cms.EDProducer('DiMuonGPUProducer@alpaka',
    muons = cms.InputTag("selectedMuons")
)

process.analyzer = cms.EDAnalyzer('DiMuonSoAAnalyzer',
                          diMuonSOA = cms.InputTag("producer"),
                          candsSOA = cms.InputTag("producer"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          HLTs = cms.vstring("HLT_DoubleMu4_3_LowMass_v"),
                          diMuonMassCut = cms.vdouble(2.5,20.0), #1S
                          keepSameSign = cms.bool(True)
)

process.dimuon_seq = cms.Sequence(
                                   process.selectedMuons *
                                   process.triggerSelection *
                                   process.producer
                                   )



process.p = cms.Path(process.dimuon_seq*process.analyzer)
