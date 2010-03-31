import FWCore.ParameterSet.Config as cms

# generator specific
from TauAnalysis.GenSimTools.gen_decaysFromZs_cfi import *

# generic
from TauAnalysis.CandidateTools.diCandidatePairProducer_cfi import diTauProducer
genZll = diTauProducer.clone()
genZll.srcLeg1 = 'genMuonsFromZs'
genZll.srcLeg2 = 'genMuonsFromZs'
genZll.srcMET = ''
genZll.recoMode = ''
genZll.useLeadingTausOnly = False

genZllReconstruction = cms.Sequence(
    produceGenDecayProductsFromZs *
    genZll
    )


