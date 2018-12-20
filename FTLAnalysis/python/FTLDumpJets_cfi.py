import FWCore.ParameterSet.Config as cms

FTLDumpJets = cms.EDAnalyzer(
    "FTLDumpJets",
    genParticlesTag = cms.untracked.InputTag("prunedGenParticles", "", ""),
    genJetsTag = cms.untracked.InputTag("slimmedGenJets", "", ""), 
    genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", ""),
    genT0Tag = cms.untracked.InputTag("genParticles", "t0", ""),
    vtxTag = cms.untracked.InputTag("offlineSlimmedPrimaryVertices", "", ""),
    vtxTimeTag = cms.untracked.InputTag("offlineSlimmedPrimaryVertices", "", ""),
    ftlRecHitsTag = cms.untracked.InputTag("ftlRecHits", "FTLBarrel", ""),
    jetsTag = cms.untracked.InputTag("slimmedJetsPuppi", "", ""),
    tracksTag = cms.untracked.InputTag("packedPFCandidates", "", ""),
    readFTLRecHits = cms.untracked.bool(True),
    treeName = cms.untracked.string("jet_tree")
)
