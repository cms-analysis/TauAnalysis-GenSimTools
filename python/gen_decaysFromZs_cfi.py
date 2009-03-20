import FWCore.ParameterSet.Config as cms

#################################################################################
#
# CV: obsolete; kept only for backwards compatibility reasons
#
decaysFromZs = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop *  ", # this is the default
    "keep+ pdgId = {Z0}",
    "drop pdgId = {Z0}"
    )
)
#################################################################################

genParticlesFromZs = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticles"),
  select = cms.vstring(
    "drop * ", # this is the default
    "keep+ pdgId = {Z0}",
    "drop pdgId = {Z0}"
  )
)

genElectronsFromZs = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromZs"),
  select = cms.vstring(
    "drop * ", # this is the default
    "keep pdgId = {e+}",
    "keep pdgId = {e-}"
  )
)

genMuonsFromZs = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromZs"),
  select = cms.vstring(
    "drop * ", # this is the default
    "keep pdgId = {mu+}",
    "keep pdgId = {mu-}"
  )
)

genTausFromZs = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromZs"),
  select = cms.vstring(
    "drop * ", # this is the default
    "keep pdgId = {tau-}",
    "keep pdgId = {tau+}"
  )
)

produceGenDecayProductsFromZs = cms.Sequence(genParticlesFromZs * genElectronsFromZs * genMuonsFromZs * genTausFromZs)


