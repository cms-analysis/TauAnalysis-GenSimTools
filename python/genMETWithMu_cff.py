import FWCore.ParameterSet.Config as cms

import RecoMET.Configuration.GenMETParticles_cff
genCandidatesForMETWithMu = RecoMET.Configuration.GenMETParticles_cff.genCandidatesForMET.clone()

genCandidatesForMETWithMu.ignoreParticleIDs = cms.vuint32(
    1000022, 2000012, 2000014,
    2000016, 1000039, 5000039,
    4000012, 9900012, 9900014,
    9900016, 39, 18, 12, 14, 16)  ###These ID's will contribute to MET

import RecoMET.METProducers.genMet_cfi
genMETWithMu = RecoMET.METProducers.genMet_cfi.genMet.clone() 
genMETWithMu.src = cms.InputTag("genCandidatesForMETWithMu") ## Output MET type
genMETWithMu.alias = cms.string('GenMETWithMu') ## Alias  for FWLite
