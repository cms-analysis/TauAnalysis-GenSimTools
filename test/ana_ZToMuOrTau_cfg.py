

import FWCore.ParameterSet.Config as cms


process = cms.Process("ANA")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )



process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        'rfio:/castor/cern.ch/user/l/lusito/SkimJanuary09/ZtautauSkimMuPFCaloTauComplStatis/muTauSkim_2.root'
        )
#                            fileNames = cms.untracked.vstring('file:PATLayer1_fromAOD_ZtoMus.root')
)

#process.load("TauAnalysis.GenSimTools.Sources.gensource_ZtoMus_PYTHIA_cfi")

# analysis
#process.load("TauAnalysis.GenSimTools.gen_cff")


# path  -----------------------------------------------------------------

#process.particleListDrawer.maxEventsToPrint = cms.untracked.int32(10)
#process.particleListDrawer.printOnlyHardInteraction = True

#process.load("TauAnalysis.GenSimTools.genDiTauReconstruction_cff")
#process.load("TauAnalysis.GenSimTools.genZllReconstruction_cff")
process.load("TauAnalysis.GenSimTools.gen_LeptondecaysFromZs_cfi")
process.include( "SimGeneral/HepPDTESSource/data/pythiapdt.cfi" )

#process.genZll.verbosity = 1



#process.printPrunedTree = cms.EDAnalyzer(
#    "ParticleTreeDrawer",
#    src = cms.InputTag("prunedGenParticles"),
#    printIndex = cms.untracked.bool(True),
#    printStatus = cms.untracked.bool(True)
#)


process.p1 = cms.Path(
#    process.genParticles + 
#    process.genParticlesPrint +
#    process.genDiTauReconstruction + 
    #process.genZllReconstruction
    process.decaysFromZs
    +process.TauMuFromZs
    +process.TauMuFromZsFilter
    #*process.printPrunedTree
)


process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *'),
    fileName = cms.untracked.string('analysis_Z.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring("p1")
        )
)

#process.load("TauAnalysis.GenSimTools.genEventContent_cff")
#process.out.outputCommands.extend( process.TauAnalysisGenSimToolsEC.outputCommands )


process.outpath = cms.EndPath( process.out )


#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('PATLayer0Summary')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    default          = cms.untracked.PSet( limit = cms.untracked.int32(0) ),
#    PATLayer0Summary = cms.untracked.PSet( limit = cms.untracked.int32(10) )
#)
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )



