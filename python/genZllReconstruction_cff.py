import FWCore.ParameterSet.Config as cms

# generator specific
from TauAnalysis.GenSimTools.gen_decaysFromZs_cfi import decaysFromZs

# generic
from TauAnalysis.CandidateTools.diTauProducer_cfi import diTauProducer
genZll = diTauProducer.clone()
genZll.hadronicTaus = 'decaysFromZs'
genZll.leptonicTaus = 'decaysFromZs'
genZll.METs = ''
genZll.metMode = 1
genZll.useLeadingTaus = False

genZllReconstruction = cms.Sequence(
    decaysFromZs *
    genZll
    )


