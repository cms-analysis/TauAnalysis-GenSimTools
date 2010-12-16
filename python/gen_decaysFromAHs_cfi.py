import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# select collections of generated electrons, muons and tau leptons
# resulting from A/H Higgs decays
#
# NOTE:
#       (1) names of particles are defined in SimGeneral/HepPDTESSource/data/pythiaparticle.tbl
#       (2) for more information about the GenParticlePruner module, see
#             https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideGenParticlePruner
#
#--------------------------------------------------------------------------------

genParticlesFromAHs = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop * ", # this is the default
        "keep+ pdgId = {h0}",                                        
        "drop pdgId = {h0}",
        "keep+ pdgId = {H0}",                                        
        "drop pdgId = {H0}",
        "keep+ pdgId = {A0}",                                        
        "drop pdgId = {A0}"                                    
    )
)

genElectronsFromAHs = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromAHs"),
  select = cms.vstring(
    "drop * ",
    "keep pdgId = {e+}",
    "keep pdgId = {e-}"
  )
)

genMuonsFromAHs = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromAHs"),
  select = cms.vstring(
    "drop * ",
    "keep pdgId = {mu+}",
    "keep pdgId = {mu-}"
  )
)

genTausFromAHs = cms.EDProducer("GenParticlePruner",
  src = cms.InputTag("genParticlesFromAHs"),
  select = cms.vstring(
    "drop * ",
    "keep pdgId = {tau+}",
    "keep pdgId = {tau-}"
  )
)

#--------------------------------------------------------------------------------
# match tau leptons resulting from A/H Higgs decays to generator level tau-jets
#--------------------------------------------------------------------------------

genTauJetsFromAHs = cms.EDProducer("TauGenJetMatchSelector",
  srcGenTauLeptons = cms.InputTag("genTausFromAHs"),
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

genElectronsFromAHtautauDecays = cms.EDFilter("TauGenJetDecayModeSelector",
    src = cms.InputTag("genTauJetsFromAHs"),
    select = cms.vstring(
        "electron"
    ),
    filter = cms.bool(False)                                       
)

genMuonsFromAHtautauDecays = cms.EDFilter("TauGenJetDecayModeSelector",
    src = cms.InputTag("genTauJetsFromAHs"),
    # list of individual tau decay modes to be used for 'select' configuation parameter
    # define in PhysicsTools/JetMCUtils/src/JetMCTag.cc
    select = cms.vstring(
        "muon"
    ),
    filter = cms.bool(False)                                       
)

genHadronsFromAHtautauDecays = cms.EDFilter("TauGenJetDecayModeSelector",
    src = cms.InputTag("genTauJetsFromAHs"),
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

genElectronsFromAHtautauDecaysWithinAcceptance = cms.EDFilter("GenJetSelector",
  src = cms.InputTag("genElectronsFromAHtautauDecays"),
  cut = cms.string('pt > 15. && abs(eta) < 2.1')
)

genMuonsFromAHtautauDecaysWithinAcceptance = cms.EDFilter("GenJetSelector",
  src = cms.InputTag("genMuonsFromAHtautauDecays"),
  cut = cms.string('pt > 15. && abs(eta) < 2.1')
)

genHadronsFromAHtautauDecaysWithinAcceptance = cms.EDFilter("GenJetSelector",
  src = cms.InputTag("genHadronsFromAHtautauDecays"),
  cut = cms.string('pt > 20. && abs(eta) < 2.3')
)

produceGenDecayProductsFromAHs = cms.Sequence(
    genParticlesFromAHs
   * genElectronsFromAHs * genMuonsFromAHs * genTausFromAHs
   * genTauJetsFromAHs  
   * genElectronsFromAHtautauDecays * genElectronsFromAHtautauDecaysWithinAcceptance
   * genMuonsFromAHtautauDecays * genMuonsFromAHtautauDecaysWithinAcceptance
   * genHadronsFromAHtautauDecays * genHadronsFromAHtautauDecaysWithinAcceptance
)


