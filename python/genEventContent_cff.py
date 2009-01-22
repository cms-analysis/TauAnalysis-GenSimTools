import FWCore.ParameterSet.Config as cms


TauAnalysisGenSimToolsEC = cms.PSet(
    outputCommands = cms.untracked.vstring(
    'keep *_tauGenJets_*_*',
    'keep *_genHadrTauSelector_*_*',
    'keep *_genLeptTauSelector_*_*',
    'keep *_genMETWithMu_*_*',
    'keep *_genMETFromNeutrals_*_*',
    'keep *_genDiTau_*_*',
    'keep *_decaysFromZs_*_*',
    'keep *_genZll_*_*'
    )
)

