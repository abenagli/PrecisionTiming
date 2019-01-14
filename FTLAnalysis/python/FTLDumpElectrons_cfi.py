import FWCore.ParameterSet.Config as cms

FTLDumpElectrons = cms.EDAnalyzer(
    'FTLDumpElectrons',
    genParticlesTag = cms.untracked.InputTag("genParticles"),
    simHitsBTLTag = cms.untracked.InputTag("g4SimHits:FastTimerHitsBarrel"),
    recHitsBTLTag = cms.untracked.InputTag("mtdRecHits:FTLBarrel"),
    clustersBTLTag = cms.untracked.InputTag("mtdClusters:FTLBarrel"),
    simHitsETLTag = cms.untracked.InputTag("g4SimHits:FastTimerHitsEndcap"),
    recHitsETLTag = cms.untracked.InputTag("mtdRecHits:FTLEndcap"),
    clustersETLTag = cms.untracked.InputTag("mtdClusters:FTLEndcap"),
    barrelElectronsTag = cms.untracked.InputTag("gedGsfElectrons"),
    endcapElectronsTag = cms.untracked.InputTag("cleanedEcalDrivenGsfElectronsFromMultiCl"),
    barrelElectronMVATag = cms.untracked.InputTag("hgcElectronMVAbarrel"),
    endcapElectronMVATag = cms.untracked.InputTag("hgcElectronMVAendcap"),
    genVtxTag = cms.untracked.InputTag("g4SimHits"),
    crysLayout = cms.untracked.int32(0),
    track_hit_DRMax = cms.double(100.05),
    track_hit_distMax = cms.double(99999.),
    treeName = cms.untracked.string("DumpElectrons"),
    verbosity = cms.bool(False)
    )
