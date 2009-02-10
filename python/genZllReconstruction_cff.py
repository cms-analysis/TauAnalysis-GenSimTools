import FWCore.ParameterSet.Config as cms

# generator specific
from TauAnalysis.GenSimTools.gen_decaysFromZs_cfi import decaysFromZs

# generic
from TauAnalysis.CandidateTools.diCandidatePairProducer_cfi import diTauProducer
genZll = diTauProducer.clone()
genZll.srcLeg1 = 'decaysFromZs'
genZll.srcLeg2 = 'decaysFromZs'
genZll.srcMET = ''
genZll.recoMode = ''
genZll.useLeadingTausOnly = False

genZllReconstruction = cms.Sequence(
    decaysFromZs *
    genZll
    )


