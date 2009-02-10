import FWCore.ParameterSet.Config as cms


from PhysicsTools.JetMCAlgos.TauGenJets_cfi import tauGenJets

# generator specific 
from TauAnalysis.GenSimTools.genMETWithMu_cff import *
from TauAnalysis.GenSimTools.genTauDecaySelector_cfi import genTauDecaySelector
from TauAnalysis.GenSimTools.genTauPairDecayFilter_cfi import genTauPairDecayFilter
# by default, genTauDecaySelector selects hadronic taus
genHadrTauSelector = genTauDecaySelector.clone()
genLeptTauSelector = genTauDecaySelector.clone()
genLeptTauSelector.tauDaughtersPdgId = (11,13)

# corrected genMet producer, it will be removed in the future
from TauAnalysis.CandidateTools.diCandidatePairProducer_cfi import diTauProducer
genDiTau = diTauProducer.clone()
genDiTau.srcLeg1 = 'genHadrTauSelector'
genDiTau.srcLeg2 = 'genLeptTauSelector'
genDiTau.srcMET = 'genMETWithMu'

genDiTauReconstruction = cms.Sequence(
    tauGenJets *
    genTauPairDecayFilter *
    genCandidatesForMETWithMu *
    genMETWithMu *
    genLeptTauSelector *
    genHadrTauSelector *
    genDiTau
    )


