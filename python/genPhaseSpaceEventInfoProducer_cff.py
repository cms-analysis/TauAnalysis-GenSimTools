import FWCore.ParameterSet.Config as cms

genPhaseSpaceEventInfo = cms.EDProducer("GenPhaseSpaceEventInfoProducer",
  srcGenEventScale = cms.InputTag('genEventScale'),
  srcGenParticles = cms.InputTag('genParticles')
)

produceGenPhaseSpaceEventInfo = cms.Sequence(genPhaseSpaceEventInfo)
