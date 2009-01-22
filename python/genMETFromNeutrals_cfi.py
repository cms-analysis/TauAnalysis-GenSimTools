import FWCore.ParameterSet.Config as cms

genMETFromNeutrals = cms.EDProducer(
    "GenMETFromNeutralsProducer",
    GenParticles = cms.InputTag('genParticles'),
    # use only neutrinos from direct taus. If set, further requirements ignored
    neutrinosFromTaus = cms.untracked.bool(False),
    minPt = cms.double(0.0),
    maxEta = cms.double(5.0), # if < 0, all etas accepted
    excludeBSM = cms.untracked.bool(False),
    verbose = cms.untracked.bool(False)
    )

