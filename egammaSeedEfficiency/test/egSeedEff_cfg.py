import FWCore.ParameterSet.Config as cms

process = cms.Process('electrons')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

## ----------------- Global Tag -----------------
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_hlt_GRun', '')

##-------------------- Report and output ---------------------------   
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.TFileService=cms.Service("TFileService",
        fileName=cms.string("output.root")
)
process.options = cms.untracked.PSet(
        allowUnscheduled = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(False),
)

##-------------------- Define the source  ----------------------------
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_12_6_0_pre4/RelValZEE_14/GEN-SIM-DIGI-RAW/125X_mcRun3_2022_realistic_v4-v1/2580000/94a518d8-73e4-4cba-a97f-23f8c2c9834c.root'
    )
)

process.egammaReconstruction = cms.EDAnalyzer(
    'egSeedingEff',
     electron = cms.InputTag('hltEgammaGsfElectrons'),
     genParticles = cms.InputTag("genParticles"),     
)

# ------------------ path --------------------------
process.p = cms.Path(process.egammaReconstruction)
