import FWCore.ParameterSet.Config as cms



genTauPairLeptHadr = cms.EDFilter(
    "GenTauPairLepHadr",
    tauGenJets = cms.InputTag('tauGenJets'),
    tauLeptonicDaughtersPdgId = cms.vint32( 13 ),
    verbose =  cms.untracked.bool(False)
    )
