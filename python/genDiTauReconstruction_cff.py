import FWCore.ParameterSet.Config as cms

## generator specific 
# genMET config (it will be removed in the future)
from TauAnalysis.GenSimTools.genMETWithMu_cff import *

# GenTauJets producer
from PhysicsTools.JetMCAlgos.TauGenJets_cfi import tauGenJets

# import config for selection of generated tau-jets
from TauAnalysis.GenSimTools.tauGenJetSelector_cfi import *

# import config for count filter
from TauAnalysis.CandidateTools.candidateCountFilter_cfi import candidateCountFilter
atLeastOneGenTauToMuon = candidateCountFilter.clone()
atLeastOneGenTauToMuon.src = 'selectedGenTauDecaysToMuon'
atLeastOneGenTauToMuon.minNumber = 1
atLeastOneGenTauToHadrons = candidateCountFilter.clone()
atLeastOneGenTauToHadrons.src = 'selectedGenTauDecaysToHadrons'
atLeastOneGenTauToHadrons.minNumber = 1

## diTau producer
from TauAnalysis.CandidateTools.diCandidatePairProducer_cfi import diTauProducer
genDiTau = diTauProducer.clone()
genDiTau.srcLeg1 = 'selectedGenTauDecaysToMuon'
genDiTau.srcLeg2 = 'selectedGenTauDecaysToHadrons'
genDiTau.srcMET = 'genMETWithMu'

genDiTauReconstruction = cms.Sequence(
    tauGenJets *
    selectedGenTauDecaysToMuon *
    selectedGenTauDecaysToHadrons *
    atLeastOneGenTauToMuon *
    atLeastOneGenTauToHadrons *
    genCandidatesForMETWithMu *
    genMETWithMu *
    genDiTau
    )
