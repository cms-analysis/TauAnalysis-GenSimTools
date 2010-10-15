import FWCore.ParameterSet.Config as cms

## generator specific 
# GenTauJets producer
from PhysicsTools.JetMCAlgos.TauGenJets_cfi import tauGenJets

# import config for selection of generated tau-jets
from TauAnalysis.GenSimTools.tauGenJetSelector_cfi import *

# import config for count filter
from TauAnalysis.CandidateTools.candidateCountFilter_cfi import candidateCountFilter
atLeastOneGenTauToElectron = candidateCountFilter.clone()
atLeastOneGenTauToElectron.src = 'selectedGenTauDecaysToElectron'
atLeastOneGenTauToElectron.minNumber = 1
atLeastOneGenTauToHadrons = candidateCountFilter.clone()
atLeastOneGenTauToHadrons.src = 'selectedGenTauDecaysToHadrons'
atLeastOneGenTauToHadrons.minNumber = 1

## diTau producer
from TauAnalysis.CandidateTools.diCandidatePairProducer_cfi import diTauProducer
genDiTau = diTauProducer.clone()
genDiTau.srcLeg1 = 'selectedGenTauDecaysToElectron'
genDiTau.srcLeg2 = 'selectedGenTauDecaysToHadrons'
genDiTau.srcMET = 'genMetTrue'

genDiTauReconstruction = cms.Sequence(
    tauGenJets *
    selectedGenTauDecaysToElectron *
    selectedGenTauDecaysToHadrons *
    atLeastOneGenTauToElectron *
    atLeastOneGenTauToHadrons *
    genDiTau
    )
