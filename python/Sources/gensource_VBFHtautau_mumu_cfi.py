import FWCore.ParameterSet.Config as cms

from Configuration.GenProduction.PythiaUESettings_cfi import *
# for every process including tau should be use TAUOLA
from TauAnalysis.GenSim.Sources.TauolaSettings_mumu_cff import *

source = cms.Source(
    "PythiaSource",
    ExternalGenerators = cms.PSet(
          Tauola = cms.untracked.PSet(
                 TauolaPolar,
                 TauolaDefaultInputCards
          ),
          parameterSets = cms.vstring('Tauola')
    ),
    UseExternalGenerators = cms.untracked.bool(True),
    PythiaParameters = cms.PSet(
       pythiaUESettingsBlock,
       myParameters = cms.vstring(
            'MSEL=0 ! Users defined processes', 
            'MSUB(123)=1             !ZZ fusion to H',
            'MSUB(124)=1             !WW fusion to H',
            'PMAS(25,1)=135.00 ! H mass', 
            'MDME(210,1)=0', 
            'MDME(211,1)=0', 
            'MDME(212,1)=0', 
            'MDME(213,1)=0', 
            'MDME(214,1)=0', 
            'MDME(215,1)=0', 
            'MDME(216,1)=0', 
            'MDME(217,1)=0', 
            'MDME(218,1)=0', 
            'MDME(219,1)=0', 
            'MDME(220,1)=1 ! switch on decay to taus', 
            'MDME(221,1)=0', 
            'MDME(222,1)=0', 
            'MDME(223,1)=0', 
            'MDME(224,1)=0', 
            'MDME(225,1)=0', 
            'MDME(226,1)=0', 
            'MSTJ(22)=2   ! Do not decay unstable particles', 
            'PARJ(71)=10. ! with c*tau > cTauMin (in mm) in PYTHIA'
            ),
       
        parameterSets = cms.vstring('pythiaUESettings','myParameters')
    )

#    maxEventsToPrint = cms.untracked.int32(-1),  
#    pythiaPylistVerbosity = cms.untracked.int32(1),  
#    pythiaHepMCVerbosity = cms.untracked.bool(False),
)


