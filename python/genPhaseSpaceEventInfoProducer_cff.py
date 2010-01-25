import FWCore.ParameterSet.Config as cms

genPhaseSpaceEventInfo = cms.EDProducer("GenPhaseSpaceEventInfoProducer",
    srcGenEventInfo = cms.InputTag('generator'),
    srcGenParticles = cms.InputTag('genParticles')
)

produceGenPhaseSpaceEventInfo = cms.Sequence(genPhaseSpaceEventInfo)
