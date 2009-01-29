import FWCore.ParameterSet.Config as cms
import copy

# import config for generator level tau-jet production
from PhysicsTools.JetMCAlgos.TauGenJets_cfi import *

# import config for selection of generated tau-jets
from TauAnalysis.GenSimTools.tauGenJetSelector_cfi import *

produceTauGenJetsForTauAnalyses = cms.Sequence( tauGenJets
                                               *selectedGenTauDecaysToElectron
                                               *selectedGenTauDecaysToElectronEta21Cumulative * selectedGenTauDecaysToElectronEta21Individual
                                               *selectedGenTauDecaysToElectronPt15Cumulative * selectedGenTauDecaysToElectronPt15Individual
                                               *selectedGenTauDecaysToMuon
                                               *selectedGenTauDecaysToMuonEta21Cumulative * selectedGenTauDecaysToMuonEta21Individual
                                               *selectedGenTauDecaysToMuonPt15Cumulative * selectedGenTauDecaysToMuonPt15Individual
                                               *selectedGenTauDecaysToHadrons
                                               *selectedGenTauDecaysToHadronsEta21Cumulative * selectedGenTauDecaysToHadronsEta21Individual
                                               *selectedGenTauDecaysToHadronsPt20Cumulative * selectedGenTauDecaysToHadronsPt20Individual )
