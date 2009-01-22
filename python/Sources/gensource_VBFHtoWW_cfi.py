import FWCore.ParameterSet.Config as cms

from Configuration.GenProduction.PythiaUESettings_cfi import *

source = cms.Source(
    "PythiaSource",
    PythiaParameters = cms.PSet(
       pythiaUESettingsBlock,
       myParameters = cms.vstring(
            'MSEL=0 ! Users defined processes', 
            'MSUB(123)=1             !ZZ fusion to H',
            'MSUB(124)=1             !WW fusion to H',
            'PMAS(25,1)=165.00 ! H mass', 
            'MDME(210,1)=0 ! H decays', 
            'MDME(211,1)=0', 
            'MDME(212,1)=0', 
            'MDME(213,1)=0', 
            'MDME(214,1)=0', 
            'MDME(215,1)=0', 
            'MDME(216,1)=0', 
            'MDME(217,1)=0', 
            'MDME(218,1)=0', 
            'MDME(219,1)=0', 
            'MDME(220,1)=0 ! H to taus', 
            'MDME(221,1)=0', 
            'MDME(222,1)=0', 
            'MDME(223,1)=0', 
            'MDME(224,1)=0', 
            'MDME(225,1)=0', 
            'MDME(226,1)=1 ! H to Ws', 
            'MDME(190,1)=0 ! W decays',
            'MDME(191,1)=0',
            'MDME(192,1)=0',
            'MDME(194,1)=0',
            'MDME(195,1)=0',
            'MDME(196,1)=0',
            'MDME(198,1)=0',
            'MDME(199,1)=0',
            'MDME(200,1)=0',
            'MDME(206,1)=1 ! W to e,nu',
            'MDME(207,1)=1 ! W to mu,nu',
            'MDME(208,1)=0 ! W to tau,nu',
            ),
       
        parameterSets = cms.vstring('pythiaUESettings','myParameters')
    )

#    maxEventsToPrint = cms.untracked.int32(-1),  
#    pythiaPylistVerbosity = cms.untracked.int32(1),  
#    pythiaHepMCVerbosity = cms.untracked.bool(False),
)
