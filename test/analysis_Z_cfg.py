

import FWCore.ParameterSet.Config as cms

# Choose analysis mode:
# 0: Z->mu,mu 
# 1: Z->tau,tau->mu,tau-jet
mode = 0

process = cms.Process("ANA")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

if mode==0:
    print "*** Z->mu,mu gen-analysis ***"
    process.load("Configuration.Generator.ZMM_cfi")
else:
    print "*** Z->tau,tau->mu,tau-jet gen-analysis ***"
    process.load("Configuration.Generator.ZTT_Tauola_OneLepton_OtherHadrons_cfi")
    process.generator.ExternalDecays.Tauola.InputCards.mdtau = 116

process.load("FastSimulation.Configuration.RandomServiceInitialization_cff")

# analysis
process.load("TauAnalysis.GenSimTools.gen_cff")


# path  -----------------------------------------------------------------

process.particleListDrawer.maxEventsToPrint = cms.untracked.int32(10)
process.particleListDrawer.printOnlyHardInteraction = True

if mode==0:
    process.load("TauAnalysis.GenSimTools.genZllReconstruction_cff")
    process.genZll.verbosity = 1
else:
    process.load("TauAnalysis.GenSimTools.genDiTauReconstruction_cff")
    process.genDiTau.verbosity = 0
    
process.p1 = cms.Path(
    process.ProductionFilterSequence + 
    process.genParticles +
    process.genParticlesPrint
    )
if mode==0:
    process.p1 += process.genZllReconstruction
else:
    process.load("RecoMET.Configuration.GenMETParticles_cff")
    process.load("RecoMET.METProducers.genMetTrue_cfi")
    process.p1 += process.genParticlesForMETAllVisible*process.genMetTrue
    process.p1 += process.genDiTauReconstruction


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

