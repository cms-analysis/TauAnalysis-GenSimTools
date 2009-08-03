import FWCore.ParameterSet.Config as cms

# The Alpgen Source. It reads unweighted alpgen files
source = cms.Source("AlpgenSource",
# use an input file name without extension unw
    fileNames = cms.untracked.vstring(
         'file:/home/cbernet/ALPGEN/v213/zjetwork/z2j'
    )
)

# The Alpgen Producer.
from GeneratorInterface.AlpgenInterface.generator_cfi import *
generator.comEnergy = 14000.0
generator.pythiaHepMCVerbosity = False
generator.maxEventsToPrint = 0
# Set the jet matching parameters as you see fit.
generator.jetMatching.applyMatching = True
generator.jetMatching.exclusive = True
generator.jetMatching.etMin = 20.0
generator.jetMatching.drMin = 0.5

# for every process including tau should be use TAUOLA
from GeneratorInterface.ExternalDecays.TauolaSettings_cff import *
generator.ExternalDecays = cms.PSet(
    Tauola = cms.untracked.PSet(
        TauolaPolar,
        InputCards = cms.PSet(
           pjak1 = cms.int32(0),
           pjak2 = cms.int32(0),
           #mdtau = cms.int32(116)  #mdtau = 0 all decays
           mdtau = cms.int32(116)  #mdtau = 116 - ONE mu+-, other taus -> all channels
        )
    ),
    parameterSets = cms.vstring('Tauola')
)

ProductionFilterSequence = cms.Sequence(generator)
