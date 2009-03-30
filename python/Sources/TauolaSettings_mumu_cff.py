import FWCore.ParameterSet.Config as cms

TauolaDefaultInputCards = cms.PSet(
    # decay all taus from hard process to muons
    # for more info see https://twiki.cern.ch/twiki/bin/view/CMS/Tauola
    InputCards = cms.vstring('TAUOLA = 0 0 102   ! TAUOLA ')
)
TauolaNoPolar = cms.PSet(
    UseTauolaPolarization = cms.bool(False)
)
TauolaPolar = cms.PSet(
    UseTauolaPolarization = cms.bool(True)
)


