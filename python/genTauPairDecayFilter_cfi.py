import FWCore.ParameterSet.Config as cms

genTauPairDecayFilter = cms.EDFilter(
    "GenTauPairDecayFilter",
    tauGenJets = cms.InputTag('tauGenJets'),
    # set pdgId of decay product or -1 for all hadronic decays
    tau1DaughtersPdgId = cms.vint32( 11, 13 ), 
    tau2DaughtersPdgId = cms.vint32( -1 ),
    verbose =  cms.untracked.bool(False)
)
