import FWCore.ParameterSet.Config as cms

FTLDumpHits = cms.EDAnalyzer(
    "FTLDumpHits",
    genParticlesTag = cms.untracked.InputTag("genParticles", "", ""),
    simHitsTag = cms.untracked.InputTag("g4SimHits", "FastTimerHitsBarrel", ""),
    tracksTag = cms.untracked.InputTag("generalTracks", "", ""),
    recHitsTag = cms.untracked.InputTag("mtdRecHits", "FTLBarrel", ""),
    treeName = cms.untracked.string("hits_tree"),
    crysLayout = cms.untracked.int32(1), # 1: tile   2: barphi   3: barz
    track_hit_DRMax = cms.double(0.05),
    track_hit_distMax = cms.double(5.00),
    verbosity = cms.bool(False)
)
