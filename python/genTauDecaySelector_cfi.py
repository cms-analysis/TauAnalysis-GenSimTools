import FWCore.ParameterSet.Config as cms

genTauDecaySelector = cms.EDFilter(
    "GenTauDecaySelector",
    src = cms.InputTag('tauGenJets'),
    # set pdgId of decay product or -1 for all hadronic decays
    tauDaughtersPdgId = cms.vint32( -1 ),
    verbose =  cms.untracked.bool(False)
)

