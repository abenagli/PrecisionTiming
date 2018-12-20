import FWCore.ParameterSet.Config as cms

FTLDumpElectronsRECO = cms.EDAnalyzer(
    "FTLDumpElectronsRECO",
    genParticlesTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    electronsTag = cms.untracked.InputTag("gedGsfElectrons", "", "HLT"),
    ftlRecHitsTag = cms.untracked.InputTag("ftlRecHits", "FTLBarrel", "HLT"),
    simTkTag = cms.untracked.InputTag("g4SimHits", "", "HLT"),
    simVtxTag = cms.untracked.InputTag("g4SimHits", "", "HLT"),
    mcTruthEleEtThr = cms.untracked.double(10),
    readFTLRecHits = cms.untracked.bool(True),
    treeName = cms.untracked.string("ele_tree")
)

FTLDumpElectronsPAT = cms.EDAnalyzer(
    "FTLDumpElectronsPAT",
    electronsTag = cms.untracked.InputTag("slimmedElectrons", "", "RECO"),
    simTkTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    simVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    mcTruthEleEtThr = cms.untracked.double(10),
    readFTLRecHits = cms.untracked.bool(True),
    treeName = cms.untracked.string("ele_tree")
)
