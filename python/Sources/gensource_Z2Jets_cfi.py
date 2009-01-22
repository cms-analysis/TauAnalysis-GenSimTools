import FWCore.ParameterSet.Config as cms

from Configuration.GenProduction.PythiaUESettings_cfi import *
# for every process including tau should be use TAUOLA
from GeneratorInterface.Pythia6Interface.TauolaSettings_cff import *


source = cms.Source("AlpgenSource",
# use input file without extension unw
    fileNames = cms.untracked.vstring(
         'file:/home/cbernet/ALPGEN/v213/zjetwork/z2j'
    ),
# for every process including tau should be use TAUOLA
    ExternalGenerators = cms.PSet(
          Tauola = cms.untracked.PSet(
                 TauolaPolar,
                 TauolaDefaultInputCards
          ),
          parameterSets = cms.vstring('Tauola')
    ),
    UseExternalGenerators = cms.untracked.bool(True),
# Pythia settings
    PythiaParameters = cms.PSet(
         pythiaUESettingsBlock,
         pythia = cms.vstring(
           'MSEL=0 ! User defined processes/Full user control',
           'MSTJ(1)=1 !...Fragmentation/hadronization on or off',
           'MSTJ(11)=3 ! Choice of the fragmentation function',
           'MSTP(143)=1 ! Call the matching routine in ALPGEN'
         ),
         parameterSets = cms.vstring( 'pythiaUESettings','pythia')
    ),
# selection of AlpGen events
    GeneratorParameters = cms.PSet(
        parameterSets = cms.vstring("generator"),
        generator = cms.vstring(
           'IXpar(2) = 1 ! inclus./exclus. sample: 0/1',
           #Inputs for clustering: minET(CLUS), deltaR(CLUS)
           'RXpar(1) = 20. ! ETCLUS : minET(CLUS)',
           'RXpar(2) = 0.7 ! RCLUS : deltaR(CLUS)'
        )
    )
)
