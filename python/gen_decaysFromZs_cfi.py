import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# select collections of generated electrons, muons and tau leptons
# resulting from Z boson decays
#
# NOTE:
#       (1) names of particles are defined in SimGeneral/HepPDTESSource/data/pythiaparticle.tbl
#       (2) for more information about the GenParticlePruner module, see
#             https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideGenParticlePruner
#
#--------------------------------------------------------------------------------

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

#--------------------------------------------------------------------------------
# match tau leptons resulting from Z boson decays to generator level tau-jets
#--------------------------------------------------------------------------------

genTauJetsFromZs = cms.EDProducer("TauGenJetMatchSelector",
  srcGenTauLeptons = cms.InputTag("genTausFromZs"),
  srcGenParticles = cms.InputTag("genParticles"),
  srcTauGenJets = cms.InputTag("tauGenJets"),
  dRmatchGenParticle = cms.double(0.1),
  dRmatchTauGenJet = cms.double(0.3)
)

#--------------------------------------------------------------------------------
# separate tau lepton decays into categories
#  o tau --> electron nu nu
#  o tau --> muon nu nu
#  o tau --> hadrons (mostly charged and neutral pions) nu
# 
# NOTE: list of individual tau decay modes to be used
#       for 'select' configuation parameter of TauGenJetDecayModeSelector module
#       define in PhysicsTools/JetMCUtils/src/JetMCTag.cc
#--------------------------------------------------------------------------------

genElectronsFromZtautauDecays = cms.EDFilter("TauGenJetDecayModeSelector",
    src = cms.InputTag("genTauJetsFromZs"),
    select = cms.vstring(
        "electron"
    ),
    filter = cms.bool(False)                                       
)

genMuonsFromZtautauDecays = cms.EDFilter("TauGenJetDecayModeSelector",
    src = cms.InputTag("genTauJetsFromZs"),
    # list of individual tau decay modes to be used for 'select' configuation parameter
    # define in PhysicsTools/JetMCUtils/src/JetMCTag.cc
    select = cms.vstring(
        "muon"
    ),
    filter = cms.bool(False)                                       
)

genHadronsFromZtautauDecays = cms.EDFilter("TauGenJetDecayModeSelector",
    src = cms.InputTag("genTauJetsFromZs"),
    select = cms.vstring(
        "oneProng0Pi0",
        "oneProng1Pi0",
        "oneProng2Pi0",
        "oneProngOther",
        "threeProng0Pi0",
        "threeProng1Pi0",
        "threeProngOther",
        "rare"
    ),
    filter = cms.bool(False)                                       
)

#--------------------------------------------------------------------------------
# select tau lepton decay products that are within detector acceptance
# (in terms of Pt and eta)
#--------------------------------------------------------------------------------

genElectronsFromZtautauDecaysWithinAcceptance = cms.EDFilter("GenJetSelector",
  src = cms.InputTag("genElectronsFromZtautauDecays"),
  cut = cms.string('pt > 15. && abs(eta) < 2.1')
)

genMuonsFromZtautauDecaysWithinAcceptance = cms.EDFilter("GenJetSelector",
  src = cms.InputTag("genMuonsFromZtautauDecays"),
  cut = cms.string('pt > 15. && abs(eta) < 2.1')
)

genHadronsFromZtautauDecaysWithinAcceptance = cms.EDFilter("GenJetSelector",
  src = cms.InputTag("genHadronsFromZtautauDecays"),
  cut = cms.string('pt > 20. && abs(eta) < 2.3')
)

produceGenDecayProductsFromZs = cms.Sequence(
    genParticlesFromZs
   * genElectronsFromZs * genMuonsFromZs * genTausFromZs
   * genTauJetsFromZs
   * genElectronsFromZtautauDecays * genElectronsFromZtautauDecaysWithinAcceptance
   * genMuonsFromZtautauDecays * genMuonsFromZtautauDecaysWithinAcceptance
   * genHadronsFromZtautauDecays * genHadronsFromZtautauDecaysWithinAcceptance
)


