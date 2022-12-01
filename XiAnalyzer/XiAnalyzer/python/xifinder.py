import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
from glob import glob

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register ('pu',
                  False, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,          # string, int, or float
                  "pu")
options.register ('new',
                  True, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,          # string, int, or float
                  "new")
options.register ('ckf',
                  True, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,          # string, int, or float
                  "ckf")

options.parseArguments()
pu=options.pu
ckf=options.ckf
new=options.new

folder = "nopu"
sample = "124X_mcRun3_2022_realistic_v7"
gt = "124X_mcRun3_2022_realistic_v7"
track="mkfit"
events = 140000

if pu:
    folder = "pu"

if ckf:
    track="ckf"

if new:
    sample = "125X_mcRun3_2022_realistic_v3"
    gt = "125X_mcRun3_2022_realistic_v3"
    events = 500000

basepath = "/lustre/cms/store/user/adiflori/xi_samples/"

path = basepath + "/" + sample + "/" + folder + "/" + track + "/*.root"

###
a0 = glob(path)
# a1 = glob('/eos/cms/store/group/phys_bphys/mkfit_pixelless_issue/xi_mkfit_nopixelless_samples/'+parts[1]+'/miniaod/mkfit/*.root')
# a2 = glob('/eos/cms/store/group/phys_bphys/mkfit_pixelless_issue/xi_mkfit_nopixelless_samples/'+parts[2]+'/miniaod/mkfit/*.root')
# a3 = glob('/eos/cms/store/group/phys_bphys/mkfit_pixelless_issue/xi_mkfit_nopixelless_samples/'+parts[3]+'/miniaod/mkfit/*.root')
# a4 = glob('/eos/cms/store/group/phys_bphys/mkfit_pixelless_issue/xi_mkfit_nopixelless_samples/'+parts[4]+'/miniaod/mkfit/*.root')
aa = ['file:'+ x for x in a0]

b0 = glob('/lustre/cms/store/user/adiflori/xi_samples/125X_mcRun3_2022_realistic_v3/nopu/ckf/*.root')
# b1 = glob('/eos/cms/store/group/phys_bphys/mkfit_pixelless_issue/xi_mkfit_nopixelless_samples/'+parts[1]+'/miniaod/mkfit_nopix/*.root')
# b2 = glob('/eos/cms/store/group/phys_bphys/mkfit_pixelless_issue/xi_mkfit_nopixelless_samples/'+parts[2]+'/miniaod/mkfit_nopix/*.root')
# b3 = glob('/eos/cms/store/group/phys_bphys/mkfit_pixelless_issue/xi_mkfit_nopixelless_samples/'+parts[3]+'/miniaod/mkfit_nopix/*.root')
# b4 = glob('/eos/cms/store/group/phys_bphys/mkfit_pixelless_issue/xi_mkfit_nopixelless_samples/'+parts[4]+'/miniaod/mkfit_nopix/*.root')
bb = ['file://'+ x for x in b0]
readFiles = cms.untracked.vstring()
readFiles.extend(aa);

print ('working on ', len(readFiles), ' files')

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(800) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(events) )

OUTNAME = sample + "_" + folder + "_" + track

print(OUTNAME)

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

## dataset status=* dataset=/ParkingBPH*/Run2018*05May2019*/AOD
process.source = cms.Source("PoolSource",
    fileNames = readFiles)
# process.source.lumisToProcess = LumiList.LumiList(filename = 'python/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_MuonPhys.txt').getVLuminosityBlockRange()
##
### JSON FILE
## https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_MuonPhys.txt


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '125X_mcRun3_2022_realistic_v3', '')## for 2016-2018ABC rereco

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.demo = cms.EDAnalyzer('XiAnalyzer',
     HLTriggerResults   = cms.InputTag("TriggerResults","","HLT"),
     beamSpotTag        = cms.InputTag("offlineBeamSpot"),
     VtxSample          = cms.InputTag('offlinePrimaryVerticesWithBS'),
     muons              = cms.InputTag('muons'), #oniaPATMuonsWithoutTrigger
     ## revtxtrks          = cms.InputTag('generalTracks'),
     fileName            = cms.untracked.string(OUTNAME+'.root'),
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string(OUTNAME+'.root')
)


process.mySequence = cms.Sequence(
                   process.demo
)

process.options.numberOfStreams = cms.untracked.uint32(0)
process.options.numberOfThreads = cms.untracked.uint32(16)

# process.dump=cms.EDAnalyzer('EventContentAnalyzer')
# process.p = cms.Path( process.dump * process.mySequence)
process.p = cms.Path( process.mySequence)
process.schedule = cms.Schedule(process.p)

'''
##########################
####### LOOPERS #########
###########################

Prompt-topology
Default: /RelValPsi2SToJPsiPiPi/CMSSW_12_0_0_pre6-120X_mcRun3_2021_realistic_v4_looper_enabled-v1/MINIAODSIM
Loopers disabled: /RelValPsi2SToJPsiPiPi/CMSSW_12_0_0_pre6-120X_mcRun3_2021_realistic_v4_looper-v1/MINIAODSIM

Displaced-topology:
Default: /RelValBuToJPsiPrimePhiToJPsiPiKK/CMSSW_12_0_0_pre6-120X_mcRun3_2021_realistic_v4_looper_enabled-v1/MINIAODSIM
Loopers disabled: /RelValBuToJPsiPrimePhiToJPsiPiKK/CMSSW_12_0_0_pre6-120X_mcRun3_2021_realistic_v4_looper-v1/MINIAODSIM


dasgoclient -query="dataset=/RelVal*/CMSSW_12_0_0_pre6-120X_mcRun3_2021_realistic_v4_looper_enabled-v1/MINIAODSIM"

###########################
###########################
###########################


###########################
#####  MKFITT  ###########
###########################

Prompt-topology:
Default: /RelValPsi2SToJPsiPiPi_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v3/MINIAODSIM
mkFit: /RelValPsi2SToJPsiPiPi_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM

DEFAULT
/RelValBuToJPsiPrimePhiToJPsiPiKK_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v2/MINIAODSIM
/RelValDisplacedSUSY_14TeV/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM
/RelValH125GGgluonfusion_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM
/RelValPsi2SToJPsiPiPi_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v3/MINIAODSIM
/RelValQCD_FlatPt_15_3000HS_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM
/RelValQCD_Pt_1800_2400_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM
/RelValSingleMuPt10/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM
/RelValSingleMuPt1000/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM
/RelValTTbar_14TeV/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM
/RelValTenTau_15_500/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM
/RelValUpsilon1SToMuMu_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM
/RelValZEE_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM
/RelValZMM_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM
/RelValZpTT_1500_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM


MKFIT:
/RelValBuToJPsiPrimePhiToJPsiPiKK_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM
/RelValDisplacedSUSY_14TeV/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM
/RelValH125GGgluonfusion_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM
/RelValPsi2SToJPsiPiPi_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM
/RelValQCD_FlatPt_15_3000HS_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM
/RelValQCD_Pt_1800_2400_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM
/RelValSingleMuPt10/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM
/RelValSingleMuPt1000/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM
/RelValTTbar_14TeV/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM
/RelValTenTau_15_500/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM
/RelValUpsilon1SToMuMu_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM
/RelValZEE_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM
/RelValZMM_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM
/RelValZpTT_1500_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM


DEFAULT ="""
dasgoclient -query="file dataset=/RelValBuToJPsiPrimePhiToJPsiPiKK_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v2/MINIAODSIM"
dasgoclient -query="file dataset=/RelValDisplacedSUSY_14TeV/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValH125GGgluonfusion_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValPsi2SToJPsiPiPi_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v3/MINIAODSIM"
dasgoclient -query="file dataset=/RelValQCD_FlatPt_15_3000HS_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValQCD_Pt_1800_2400_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValSingleMuPt10/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValSingleMuPt1000/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValTTbar_14TeV/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValTenTau_15_500/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValUpsilon1SToMuMu_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValZEE_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValZMM_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValZpTT_1500_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1-v1/MINIAODSIM"
"""
'/store/relval/CMSSW_12_1_0_pre2/RelValBuToJPsiPrimePhiToJPsiPiKK_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v2/10000/cf6880cf-5edb-49ff-9755-4ca06b4ff12e.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValDisplacedSUSY_14TeV/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v1/10000/bbaf40e4-a7a9-463e-9ab9-ce88be006473.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValH125GGgluonfusion_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v1/10000/cb80d5ff-daf5-492c-ac45-8d6127c4ae5b.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValPsi2SToJPsiPiPi_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v3/10000/934b45ff-3f99-4f69-9648-bfd129ab7d90.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v1/10000/1c6ace13-ca4d-4441-99ec-25ff4f8a9a21.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValQCD_Pt_1800_2400_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v1/10000/3315d8ce-fbe0-4d14-865a-b226ae926732.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValSingleMuPt10/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v1/10000/0c35330c-275c-48a4-8742-059f3cc2f310.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValSingleMuPt1000/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v1/10000/c3fffcdc-c393-4c10-bdfc-a98e4eda2f7e.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValTTbar_14TeV/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v1/10000/78d04134-4010-4c4d-9ed2-6cf597f0f554.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValTenTau_15_500/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v1/10000/5969c449-8b8e-4217-9295-bec361253759.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValUpsilon1SToMuMu_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v1/10000/6215ebd0-482c-46f1-99e8-6a56275e6e36.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValZEE_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v1/10000/eda11efa-97e5-4be5-9e33-345eaaf1f0cc.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValZMM_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v1/10000/67462d71-289e-4e98-a133-c607a5da3ebb.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValZpTT_1500_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v1/10000/c2d8907c-0ad1-492a-bf09-7552e206c0b5.root',



MKFIT = """
dasgoclient -query="file dataset=/RelValBuToJPsiPrimePhiToJPsiPiKK_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValDisplacedSUSY_14TeV/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValH125GGgluonfusion_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValPsi2SToJPsiPiPi_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValQCD_FlatPt_15_3000HS_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValQCD_Pt_1800_2400_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValSingleMuPt10/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValSingleMuPt1000/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValTTbar_14TeV/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValTenTau_15_500/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValUpsilon1SToMuMu_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValZEE_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValZMM_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
dasgoclient -query="file dataset=/RelValZpTT_1500_14/CMSSW_12_1_0_pre2-121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/MINIAODSIM"
"""
'/store/relval/CMSSW_12_1_0_pre2/RelValBuToJPsiPrimePhiToJPsiPiKK_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/c8a9d16d-dfcf-44a5-b651-7a86973492d2.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValDisplacedSUSY_14TeV/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/ec4aeeb4-534e-4e5d-973f-b02e12c1566f.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValH125GGgluonfusion_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/36c4035f-2fe4-4d6c-bed0-04069b28dca8.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValPsi2SToJPsiPiPi_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/3a2d63c2-d51b-40c6-b2a3-79c54b4ff4b6.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValQCD_FlatPt_15_3000HS_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/107e3906-aea5-48e9-9442-2a7cf9e94008.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValQCD_Pt_1800_2400_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/f6008e71-7e75-4e0b-80f9-e33c335c81a7.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValSingleMuPt10/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/67452eae-e7aa-49f5-a9af-956a5b60e9f1.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValSingleMuPt1000/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/6a6d80d0-73a8-4d42-9dc0-97aeeb9ff7bd.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValTTbar_14TeV/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/33c870c2-b7d8-437f-aea1-74b1e0629783.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValTenTau_15_500/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/77bc64ac-c8c0-479e-804d-1937cdaaef0d.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValUpsilon1SToMuMu_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/4d345516-d0e9-4fe3-a016-0cacde8a647a.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValZEE_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/cde5ad9e-4fba-4502-a226-38224bab2528.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValZMM_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/a9875947-527c-405e-9732-e3aa4c9d56ab.root',
'/store/relval/CMSSW_12_1_0_pre2/RelValZpTT_1500_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1_TkmkFitRecoOnly-v1/10000/fc967cf3-7373-461d-8f1d-ff144602e0fc.root',

###########################
###########################
###########################

'''
