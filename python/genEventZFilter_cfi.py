import FWCore.ParameterSet.Config as cms



genEventZFilter = cms.EDFilter(
    "GenEventZFilter",
    GenParticles = cms.InputTag('genParticles'),
    Z0DaughtersPdgId = cms.vint32(),
    verbose =  cms.untracked.bool(False)
    )
