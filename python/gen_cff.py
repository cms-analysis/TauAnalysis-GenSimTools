# event generation
# analysis of the gen event

import FWCore.ParameterSet.Config as cms


from FastSimulation.Configuration.RandomServiceInitialization_cff import *
from PhysicsTools.HepMCCandAlgos.genParticles_cfi import * 
from SimGeneral.HepPDTESSource.pythiapdt_cfi import * 
#COLIN: not allowed here 
#from TauAnalysis.GenSimTools.genEventAnalysisProducer_cfi import * 


particleListDrawer = cms.EDAnalyzer(
    "ParticleListDrawer",
    printOnlyHardInteraction = cms.untracked.bool(True),
    maxEventsToPrint = cms.untracked.int32(-1),
    src = cms.InputTag('genParticles')
  )

particleTreeDrawer = cms.EDAnalyzer(
    "ParticleTreeDrawer",
    src = cms.InputTag('genParticles'),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi =  cms.untracked.bool(True),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(False),
    printIndex = cms.untracked.bool(False),
    status = cms.untracked.vint32( 3, 2 )
  )

particleTreeDrawer3 = cms.EDAnalyzer(
    "ParticleTreeDrawer",
    src = cms.InputTag('genParticles'),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi =  cms.untracked.bool(True),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(False),
    printIndex = cms.untracked.bool(False),
    status = cms.untracked.vint32( 3 )
  )

# VBFHTauTau analysis

genParticlesPrint = cms.Sequence(
    particleListDrawer 
    )

#COlIN this is not allowed to be here. However we
# could use
#  the pruner to get the tau daughters,
#  the tau gen jet producer
#genParticlesAnalysis = cms.Sequence(
#    genParticlesPrint * 
#    genEventAnalysisProducer
#    )

genoutput = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('gen.root')
)
outpath = cms.EndPath(genoutput)


