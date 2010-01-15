import FWCore.ParameterSet.Config as cms

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
    "drop * ",
    "keep pdgId = {e+}",
    "keep pdgId = {e-}"
  )
)

genMuonsFromZs = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromZs"),
  select = cms.vstring(
    "drop * ",
    "keep pdgId = {mu+}",
    "keep pdgId = {mu-}"
  )
)

genTausFromZs = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromZs"),
  select = cms.vstring(
    "drop * ",
    "keep pdgId = {tau+}",
    "keep pdgId = {tau-}"
  )
)

genParticlesFromTauonicZdecays = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genTausFromZs"),
  select = cms.vstring(
    "drop * ",
    "keep+ pdgId = {tau+}",
    "keep+ pdgId = {tau-}",
    "drop pdgId = {tau+}",
    "drop pdgId = {tau-}"
  )
)

genElectronsFromTauonicZdecays = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromTauonicZdecays"),
  select = cms.vstring(
    "drop * ",
    "keep pdgId = {e+}",
    "keep pdgId = {e-}",
  )
)

genMuonsFromTauonicZdecays = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromTauonicZdecays"),
  select = cms.vstring(
    "drop * ",
    "keep pdgId = {mu+}",
    "keep pdgId = {mu-}"
  )
)

genHadronsFromTauonicZdecays = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromTauonicZdecays"),
  select = cms.vstring(
    "keep * ",
    "drop pdgId = {nu_tau}",
    "drop pdgId = {nu_taubar}",
    "drop pdgId = {mu+}",
    "drop pdgId = {mu-}",
    "drop pdgId = {nu_mu}",
    "drop pdgId = {nu_mubar}",
    "drop pdgId = {e+}",
    "drop pdgId = {e-}",
    "drop pdgId = {nu_e}",
    "drop pdgId = {nu_ebar}"
  )
)

produceGenDecayProductsFromZs = cms.Sequence(
    genParticlesFromZs
   * genElectronsFromZs * genMuonsFromZs * genTausFromZs
   * genParticlesFromTauonicZdecays
   * genElectronsFromTauonicZdecays * genMuonsFromTauonicZdecays * genHadronsFromTauonicZdecays
)


