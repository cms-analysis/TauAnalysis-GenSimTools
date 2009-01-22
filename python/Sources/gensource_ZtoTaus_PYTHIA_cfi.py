import FWCore.ParameterSet.Config as cms

from Configuration.GenProduction.PythiaUESettings_cfi import *
# for every process including tau should be use TAUOLA
from GeneratorInterface.Pythia6Interface.TauolaSettings_cff import *

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
    pythiaVerbosity = cms.untracked.bool(False),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        pythiaZtt = cms.vstring(
            'MSEL = 11         !Z production', 
            'MDME( 174,1) = 0  !Z decay into d dbar', 
            'MDME( 175,1) = 0  !Z decay into u ubar', 
            'MDME( 176,1) = 0  !Z decay into s sbar', 
            'MDME( 177,1) = 0  !Z decay into c cbar', 
            'MDME( 178,1) = 0  !Z decay into b bbar', 
            'MDME( 179,1) = 0  !Z decay into t tbar', 
            'MDME( 182,1) = 0  !Z decay into e- e+', 
            'MDME( 183,1) = 0  !Z decay into nu_e nu_ebar', 
            'MDME( 184,1) = 0  !Z decay into mu- mu+', 
            'MDME( 185,1) = 0  !Z decay into nu_mu nu_mubar', 
            'MDME( 186,1) = 1  !Z decay into tau- tau+', 
            'MDME( 187,1) = 0  !Z decay into nu_tau nu_taubar', 
            'MSTJ(22)=2        ! Do not decay unstable particles', 
            'PARJ(71)=10.      ! with c*tau > cTauMin (in mm) in PYTHIA'
            ),
        parameterSets = cms.vstring('pythiaUESettings','pythiaZtt')
    )
)


