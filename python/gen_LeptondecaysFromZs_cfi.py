import FWCore.ParameterSet.Config as cms

###my idea is to select first Zs that decay into leptons  --->producer+filter

decaysFromZs = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop *  ", # this is the default
    "keep+ pdgId = {Z0}",
    "drop pdgId = {Z0}"
    )
)

DecaysFromZsFilter = cms.EDFilter(
    "PATCandViewCountFilter",#PATCandViewCountFilter PATObjectFilter
     src = cms.InputTag("decaysFromZs"),
     minNumber = cms.uint32(0),
     maxNumber = cms.uint32(9999),
     #filter = cms.bool(True)
)


### then select the tau's and muon's from the decay products

TauMuFromZs = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop *  ", # this is the default
    "keep+ pdgId = {mu+}",
    "keep+ pdgId = {mu-}",
    "keep+ pdgId = {tau-}",
    "keep+ pdgId = {tau+}",
    
    )
)

TauMuFromZsFilter = cms.EDFilter(
    "PATCandViewCountFilter",#
     src = cms.InputTag("TauMuFromZs"),
     minNumber = cms.uint32(0),
     maxNumber = cms.uint32(9999),
     #filter = cms.bool(True)
)

### and then put the sequence of this filters in negation in order to reject events that contains Z's which decays into muons or taus
