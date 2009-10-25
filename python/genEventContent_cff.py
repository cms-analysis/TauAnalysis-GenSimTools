import FWCore.ParameterSet.Config as cms


TauAnalysisGenSimToolsEC = cms.PSet(
    outputCommands = cms.untracked.vstring(
    'keep *_tauGenJets_*_*',
    'keep *_selectedGenTauDecaysToMuon_*_*',
    'keep *_selectedGenTauDecaysToElectron_*_*',
    'keep *_selectedGenTauDecaysToHadrons_*_*',
    'keep *_genMetTrue_*_*',
    'keep *_genMETFromNeutrals_*_*',
    'keep *_genDiTau_*_*',
    'keep *_decaysFromZs_*_*',
    'keep *_genZll_*_*'
    )
)

