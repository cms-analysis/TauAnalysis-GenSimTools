

import FWCore.ParameterSet.Config as cms


process = cms.Process("ANA")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.load("Configuration.Generator.ZMM_cfi")
process.load("FastSimulation.Configuration.RandomServiceInitialization_cff")

# analysis
process.load("TauAnalysis.GenSimTools.gen_cff")


# path  -----------------------------------------------------------------

process.particleListDrawer.maxEventsToPrint = cms.untracked.int32(10)
process.particleListDrawer.printOnlyHardInteraction = True

process.load("TauAnalysis.GenSimTools.genZllReconstruction_cff")

process.genZll.verbosity = 1

process.p1 = cms.Path(
    process.ProductionFilterSequence + 
    process.genParticles +
    process.genParticlesPrint +
    process.genZllReconstruction
    )


process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *'),
    fileName = cms.untracked.string('analysis_Z.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring("p1")
        )
)

process.load("TauAnalysis.GenSimTools.genEventContent_cff")
process.out.outputCommands.extend( process.TauAnalysisGenSimToolsEC.outputCommands )


process.outpath = cms.EndPath( process.out )


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )



