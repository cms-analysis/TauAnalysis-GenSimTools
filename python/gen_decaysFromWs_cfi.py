import FWCore.ParameterSet.Config as cms

from SimGeneral.HepPDTESSource.pythiapdt_cfi import *

#--------------------------------------------------------------------------------
# select collections of generated electrons, muons and tau leptons
# resulting from W boson decays
#
# NOTE:
#       (1) names of particles are defined in SimGeneral/HepPDTESSource/data/pythiaparticle.tbl
#       (2) for more information about the GenParticlePruner module, see
#             https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideGenParticlePruner
#
#--------------------------------------------------------------------------------

genParticlesFromWs = cms.EDProducer("GenParticlesFromZsSelector",
    src = cms.InputTag("genParticles"),
    pdgIdsMothers = cms.vint32(24),
    pdgIdsDaughters = cms.vint32(11, 12, 13, 14, 15, 16),
    maxDaughters = cms.int32(2),
    minDaughters = cms.int32(2)
)

genElectronsFromWs = cms.EDProducer("GenParticlesFromZsSelector",
    src = cms.InputTag("genParticles"),
    pdgIdsMothers = cms.vint32(24),
    pdgIdsDaughters = cms.vint32(11),
    maxDaughters = cms.int32(1),
    minDaughters = cms.int32(1)
)

genMuonsFromWs = cms.EDProducer("GenParticlesFromZsSelector",
    src = cms.InputTag("genParticles"),
    pdgIdsMothers = cms.vint32(24),
    pdgIdsDaughters = cms.vint32(13),
    maxDaughters = cms.int32(1),
    minDaughters = cms.int32(1)
)

genTausFromWs = cms.EDProducer("GenParticlesFromZsSelector",
    src = cms.InputTag("genParticles"),
    pdgIdsMothers = cms.vint32(24),
    pdgIdsDaughters = cms.vint32(15),
    maxDaughters = cms.int32(1),
    minDaughters = cms.int32(1)
)

genNusFromWs = cms.EDProducer("GenParticlesFromZsSelector",
    src = cms.InputTag("genParticles"),
    pdgIdsMothers = cms.vint32(24),
    pdgIdsDaughters = cms.vint32(12, 14, 16),
    maxDaughters = cms.int32(2),
    minDaughters = cms.int32(1)
)

#--------------------------------------------------------------------------------
# match tau leptons resulting from W boson decays to generator level tau-jets
#--------------------------------------------------------------------------------

genTauJetsFromWs = cms.EDProducer("TauGenJetMatchSelector",
    srcGenTauLeptons = cms.InputTag("genTausFromWs"),
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
#       defined in PhysicsTools/JetMCUtils/src/JetMCTag.cc
#--------------------------------------------------------------------------------

genElectronsFromWtaunuDecays = cms.EDFilter("TauGenJetDecayModeSelector",
    src = cms.InputTag("genTauJetsFromWs"),
    select = cms.vstring("electron"),
    filter = cms.bool(False)
)

genMuonsFromWtaunuDecays = cms.EDFilter("TauGenJetDecayModeSelector",
    src = cms.InputTag("genTauJetsFromWs"),
    select = cms.vstring("muon"),
    filter = cms.bool(False)
)

genHadronsFromWtaunuDecays = cms.EDFilter("TauGenJetDecayModeSelector",
    src = cms.InputTag("genTauJetsFromWs"),
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
# select decay products that are within detector acceptance (in terms of Pt and eta)
#--------------------------------------------------------------------------------

genHadronsFromWtaunuDecaysWithinAcceptance = cms.EDFilter("GenJetSelector",
    src = cms.InputTag("genHadronsFromWtaunuDecays"),
    cut = cms.string('pt > 30 & abs(eta) < 2.3')
)

genNusFromWsWithinAcceptance = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genNusFromWs"),
    select = cms.vstring('keep pt > 35')
)
                                                  
produceGenDecayProductsFromWs = cms.Sequence(
    genParticlesFromWs    
   * genElectronsFromWs * genMuonsFromWs * genTausFromWs * genNusFromWs
   * genTauJetsFromWs
   * genElectronsFromWtaunuDecays * genMuonsFromWtaunuDecays * genHadronsFromWtaunuDecays
   * genHadronsFromWtaunuDecaysWithinAcceptance * genNusFromWsWithinAcceptance
)
