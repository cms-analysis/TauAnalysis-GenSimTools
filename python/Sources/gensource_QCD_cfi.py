import FWCore.ParameterSet.Config as cms

from Configuration.GenProduction.PythiaUESettings_cfi import *

source = cms.Source(
    "PythiaSource",
    PythiaParameters = cms.PSet(
       pythiaUESettingsBlock,
       myParameters = cms.vstring(
         'MSEL=1               ! QCD hight pT processes',
         'CKIN(3)=50.          ! minimum pt hat for hard interactions',
         'CKIN(4)=80.          ! maximum pt hat for hard interactions'    
       ),
       parameterSets = cms.vstring('pythiaUESettings','myParameters')
    )
)


