import FWCore.ParameterSet.Config as cms

TauolaDefaultInputCards = cms.PSet(
    # decay one tau from hard process to muon other taus to hadrons
    # for more info see https://twiki.cern.ch/twiki/bin/view/CMS/Tauola
    InputCards = cms.vstring('TAUOLA = 0 0 116   ! TAUOLA ')
)
TauolaNoPolar = cms.PSet(
    UseTauolaPolarization = cms.bool(False)
)
TauolaPolar = cms.PSet(
    UseTauolaPolarization = cms.bool(True)
)


