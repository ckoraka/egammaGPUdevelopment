import FWCore.ParameterSet.Config as cms

Ntuplizer = cms.EDAnalyzer("Ntuplizer",
    treeName = cms.string("TagAndProbe"),
    electrons = cms.InputTag("gedGsfElectrons"),
    genParticles = cms.InputTag("genParticles"), 
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp80"),
    eleLooseIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wpLoose"),
    triggerSet     = cms.InputTag("patTriggerUnpacker"),
    triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
    L1EG = cms.InputTag("rawCaloStage2Digis","EGamma"),
    L1EmuEG = cms.InputTag("simCaloStage2Digis",""),
    recHitsEB = cms.InputTag("ecalRecHit:EcalRecHitsEB"),
    recHitsEE = cms.InputTag("ecalRecHit:EcalRecHitsEE"),
    calotower = cms.InputTag("rawCaloStage2Digis","CaloTower"),
    uncalibratedRecHitCollectionEB = cms.InputTag("ecalMaxSampleUncalibRecHit","EcalUncalibRecHitsEB"),
    uncalibratedRecHitCollectionEE = cms.InputTag("ecalMaxSampleUncalibRecHit","EcalUncalibRecHitsEE"),
    GTRecordCollection = cms.string('gtDigis'),
    Vertices = cms.InputTag("offlinePrimaryVertices"),
    triggerListTag = HLTLISTTAG,
    triggerListProbe = HLTLISTPROBE,
    useGenMatch = cms.bool(False),
    useHLTMatch = cms.bool(False)
)

NtupleSeq = cms.Sequence(
    hltFilter           +
    patTriggerUnpacker  +
    Ntuplizer
)

