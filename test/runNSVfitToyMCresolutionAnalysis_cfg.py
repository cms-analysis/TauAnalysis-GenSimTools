import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('runNSVfitToyMCresolutionAnalysis')

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.MessageLogger.suppressInfo = cms.untracked.vstring()
process.MessageLogger.suppressWarning = cms.untracked.vstring("PATTriggerProducer",)
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START311_V2::All')

#--------------------------------------------------------------------------------
# produce "toy" Monte Carlo pat::Electron, pat::Muon and pat::Tau objects

process.load('PhysicsTools/JetMCAlgos/TauGenJets_cfi')
process.load('TauAnalysis/GenSimTools/gen_decaysFromZs_cfi')
#process.load('TauAnalysis/GenSimTools/gen_decaysFromAHs_cfi.py')

process.convertGenElectrons = cms.EDProducer("GenParticleFromTauGenJetProducer",
    src = cms.InputTag('genElectronsFromZtautauDecays')
)

process.produceToyElectrons = cms.EDProducer("ToyPATElectronProducer",
    src = cms.InputTag('convertGenElectrons'),
    resolution = cms.double(0.)
)

process.selectToyElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag('produceToyElectrons'),                   
    cut = cms.string('pt > 10. & abs(eta) < 2.'),
    filter = cms.bool(False)
)

process.smearToyElectrons = cms.EDProducer("ToyPATElectronProducer",
    src = cms.InputTag('convertGenElectrons'),
    resolution = cms.double(0.03)
)

process.selectSmearedToyElectrons = process.selectToyElectrons.clone(
    src = cms.InputTag('smearToyElectrons'),                   
)

process.convertGenMuons = cms.EDProducer("GenParticleFromTauGenJetProducer",
    src = cms.InputTag('genMuonsFromZtautauDecays')
)

process.produceToyMuons = cms.EDProducer("ToyPATMuonProducer",
    src = cms.InputTag('convertGenMuons'),
    resolution = cms.double(0.)
)

process.selectToyMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag('produceToyMuons'),                   
    cut = cms.string('pt > 10. & abs(eta) < 2.'),
    filter = cms.bool(False)
)

process.smearToyMuons = cms.EDProducer("ToyPATMuonProducer",
    src = cms.InputTag('convertGenMuons'),
    resolution = cms.double(0.03)
)

process.selectSmearedToyMuons = process.selectToyMuons.clone(
    src = cms.InputTag('smearToyMuons'),                   
)

process.produceToyTaus = cms.EDProducer("ToyPATTauProducer",
    src = cms.InputTag('genHadronsFromZtautauDecays'),
    resolution = cms.double(0.)
)

process.selectToyTaus = cms.EDFilter("PATTauSelector",
    src = cms.InputTag('produceToyTaus'),
    cut = cms.string('pt > 10. & abs(eta) < 2.'),
    filter = cms.bool(False)
)

process.smearToyTaus = cms.EDProducer("ToyPATTauProducer",
    src = cms.InputTag('genHadronsFromZtautauDecays'),
    resolution = cms.double(0.10)
)

process.selectSmearedToyTaus = process.selectToyTaus.clone(
    src = cms.InputTag('smearToyTaus'),                   
)

process.produceToyMEt = cms.EDProducer("ToyPATMEtProducer",
    srcGenParticles = cms.InputTag('genTausFromZs'),  
    resolutionX = cms.double(0.),
    resolutionY = cms.double(0.)
)

process.produceSmearedToyMEt = cms.EDProducer("ToyPATMEtProducer",
    srcGenParticles = cms.InputTag('genTausFromZs'),
    srcElectrons = cms.InputTag('selectToyElectrons'), 
    srcMuons = cms.InputTag('selectToyMuons'), 
    srcTaus = cms.InputTag('selectToyTaus'),                                       
    resolutionX = cms.double(5.0),
    resolutionY = cms.double(5.0)
)

process.produceToyParticles = cms.Sequence(
    process.tauGenJets
   + process.produceGenDecayProductsFromZs
   + process.convertGenElectrons + process.produceToyElectrons + process.selectToyElectrons
   + process.smearToyElectrons + process.selectSmearedToyElectrons
   + process.convertGenMuons + process.produceToyMuons + process.selectToyMuons
   + process.smearToyMuons + process.selectSmearedToyMuons
   + process.produceToyTaus + process.selectToyTaus
   + process.smearToyTaus + process.selectSmearedToyTaus
   + process.produceToyMEt + process.produceSmearedToyMEt
)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# select MET > 10 GeV events

process.filterToyMEt = cms.EDFilter("PATMETSelector",
    src = cms.InputTag('produceToyMEt'),                   
    cut = cms.string('pt > 10.'),
    filter = cms.bool(True)
)

process.filterSmearedToyMEt = cms.EDFilter("PATMETSelector",
    src = cms.InputTag('produceSmearedToyMEt'),                   
    cut = cms.string('pt > 10.'),
    filter = cms.bool(True)
)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# run NSVfit algorithm for combinations of
#  o  e  +  mu pairs
#  o  mu + tau pairs
#  o tau + tau pairs

process.load('TauAnalysis/CandidateTools/nSVfitAlgorithmDiTau_cfi')
from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName

process.runIdealResolution = cms.Sequence(process.filterToyMEt)
process.runRealisticResolution = cms.Sequence(process.filterSmearedToyMEt)

nSVfitRESoptions = {
    'idealResolution' : {
        'sequence'       : process.runIdealResolution,
        'srcElectrons'   : cms.InputTag('selectToyElectrons'),
        'srcMuons'       : cms.InputTag('selectToyMuons'),
        'srcTaus'        : cms.InputTag('selectToyTaus'),
        'srcMEt'         : cms.InputTag('filterToyMEt'),
        'metResolutionX' : cms.string("0.5"),
        'metResolutionY' : cms.string("0.5")
    },
    'realisticResolution' : {
        'sequence'       : process.runRealisticResolution,
        'srcElectrons'   : cms.InputTag('selectSmearedToyElectrons'),
        'srcMuons'       : cms.InputTag('selectSmearedToyMuons'),
        'srcTaus'        : cms.InputTag('selectSmearedToyTaus'),
        'srcMEt'         : cms.InputTag('filterSmearedToyMEt'),
        'metResolutionX' : cms.string("5.0"),
        'metResolutionY' : cms.string("5.0")
    }
}

nSVfitLLoptions = {
    'psKineMEt'     : cms.VPSet(),
    'psKineMEtLogM' : cms.VPSet(process.nSVfitResonanceLikelihoodLogM)
}    

#nSVfitProducer_template = process.nSVfitProducerByIntegration
nSVfitProducer_template = process.nSVfitProducerByLikelihoodMaximization

for optRESnSVfitName, optRESnSVfit in nSVfitRESoptions.items():

    # adjust MET likelihood to "toy" Monte Carlo resolution
    nSVfitEventLikelihoodMEt_customized = cms.PSet(
        pluginName = cms.string("nSVfitEventLikelihoodMEt"),
        pluginType = cms.string("NSVfitEventLikelihoodMEt"),
        resolution = cms.PSet(
            parSigma = optRESnSVfit['metResolutionX'],
            parBias = cms.string("0."),
            perpSigma = optRESnSVfit['metResolutionY'],
            perpBias = cms.string("0.")
        )
    )
    setattr(process, composeModuleName(["nSVfitEventLikelihoodMEt_customized", optRESnSVfitName]), nSVfitEventLikelihoodMEt_customized)
        
    for optLLnSVfitName, optLLnSVfit in nSVfitLLoptions.items():
                
        nSVfitElecMuPairHypotheses = copy.deepcopy(nSVfitProducer_template)
        nSVfitElecMuPairHypotheses.config.event.resonances.A.daughters.leg1 = \
          process.nSVfitConfig_template.event.resonances.A.daughters.leg1.clone(
            src = optRESnSVfit['srcElectrons'],
            likelihoodFunctions = cms.VPSet(process.nSVfitElectronLikelihoodPhaseSpace),
            builder = process.nSVfitTauToElecBuilder
        )
        nSVfitElecMuPairHypotheses.config.event.resonances.A.daughters.leg2 = \
          process.nSVfitConfig_template.event.resonances.A.daughters.leg2.clone(
            src = optRESnSVfit['srcMuons'],
            likelihoodFunctions = cms.VPSet(process.nSVfitMuonLikelihoodPhaseSpace),
            builder = process.nSVfitTauToMuBuilder
        )
        nSVfitElecMuPairHypotheses.config.event.resonances.A.likelihoodFunctions = optLLnSVfit
        nSVfitElecMuPairHypotheses.config.event.srcMEt = optRESnSVfit['srcMEt']
        nSVfitElecMuPairHypotheses.config.event.likelihoodFunctions = cms.VPSet(nSVfitEventLikelihoodMEt_customized)
        setattr(process, composeModuleName(["nSVfitElecMuPairHypotheses", optRESnSVfitName, optLLnSVfitName]), nSVfitElecMuPairHypotheses)
        optRESnSVfit['sequence'] += nSVfitElecMuPairHypotheses

        nSVfitMuTauPairHypotheses = copy.deepcopy(nSVfitProducer_template)
        nSVfitMuTauPairHypotheses.config.event.resonances.A.daughters.leg1 = \
          process.nSVfitConfig_template.event.resonances.A.daughters.leg1.clone(
            src = optRESnSVfit['srcMuons'],
            likelihoodFunctions = cms.VPSet(process.nSVfitMuonLikelihoodPhaseSpace),
            builder = process.nSVfitTauToMuBuilder
        )
        nSVfitMuTauPairHypotheses.config.event.resonances.A.daughters.leg2 = \
          process.nSVfitConfig_template.event.resonances.A.daughters.leg2.clone(
            src = optRESnSVfit['srcTaus'],
            likelihoodFunctions = cms.VPSet(process.nSVfitTauLikelihoodPhaseSpace),
            builder = process.nSVfitTauToHadBuilder
        )
        nSVfitMuTauPairHypotheses.config.event.resonances.A.likelihoodFunctions = optLLnSVfit
        nSVfitMuTauPairHypotheses.config.event.srcMEt = optRESnSVfit['srcMEt']
        nSVfitMuTauPairHypotheses.config.event.likelihoodFunctions = cms.VPSet(nSVfitEventLikelihoodMEt_customized)
        setattr(process, composeModuleName(["nSVfitMuTauPairHypotheses", optRESnSVfitName, optLLnSVfitName]), nSVfitMuTauPairHypotheses)
        optRESnSVfit['sequence'] += nSVfitMuTauPairHypotheses

        nSVfitDiTauPairHypotheses = copy.deepcopy(nSVfitProducer_template)
        nSVfitDiTauPairHypotheses.config.event.resonances.A.daughters.leg1 = \
          process.nSVfitConfig_template.event.resonances.A.daughters.leg1.clone(
            src = optRESnSVfit['srcTaus'],
            likelihoodFunctions = cms.VPSet(process.nSVfitTauLikelihoodPhaseSpace),
            builder = process.nSVfitTauToHadBuilder
        )
        nSVfitDiTauPairHypotheses.config.event.resonances.A.daughters.leg2 = \
          process.nSVfitConfig_template.event.resonances.A.daughters.leg2.clone(
            src = optRESnSVfit['srcTaus'],
            likelihoodFunctions = cms.VPSet(process.nSVfitTauLikelihoodPhaseSpace),
            builder = process.nSVfitTauToHadBuilder
        )
        nSVfitDiTauPairHypotheses.config.event.resonances.A.likelihoodFunctions = optLLnSVfit
        nSVfitDiTauPairHypotheses.config.event.srcMEt = optRESnSVfit['srcMEt']
        nSVfitDiTauPairHypotheses.config.event.likelihoodFunctions = cms.VPSet(nSVfitEventLikelihoodMEt_customized)
        setattr(process, composeModuleName(["nSVfitDiTauPairHypotheses", optRESnSVfitName, optLLnSVfitName]), nSVfitDiTauPairHypotheses)
        optRESnSVfit['sequence'] += nSVfitDiTauPairHypotheses

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# fill hisrograms of NSVfit reconstructed mass
# and control plots for "toy" Monte Carlo pat::Electron, pat::Muon and pat::Tau objects

for optRESnSVfitName, optRESnSVfit in nSVfitRESoptions.items():
    for optLLnSVfitName, optLLnSVfit in nSVfitLLoptions.items():
            
        analyzeElecMuPairs = cms.EDAnalyzer("ToyMCnSVfitElecMuPairResolutionAnalyzer",
            srcNSVfitEventHypothesis = cms.InputTag(composeModuleName(["nSVfitElecMuPairHypotheses", optRESnSVfitName, optLLnSVfitName])),
            srcGenParticles = cms.InputTag('genTausFromZs'),
            dqmDirectory = cms.string("/%s/nSVfitElecMuPairs/%s" % (optRESnSVfitName, optLLnSVfitName))
        )
        setattr(process, composeModuleName(["analyzeElecMuPairs", optRESnSVfitName, optLLnSVfitName]), analyzeElecMuPairs)
        optRESnSVfit['sequence'] += analyzeElecMuPairs
    
        analyzeMuTauPairs = cms.EDAnalyzer("ToyMCnSVfitMuTauPairResolutionAnalyzer",
            srcNSVfitEventHypothesis = cms.InputTag(composeModuleName(["nSVfitMuTauPairHypotheses", optRESnSVfitName, optLLnSVfitName])),
            srcGenParticles = cms.InputTag('genTausFromZs'),
            dqmDirectory = cms.string("/%s/nSVfitMuTauPairs/%s" % (optRESnSVfitName, optLLnSVfitName))
        )
        setattr(process, composeModuleName(["analyzeMuTauPairs", optRESnSVfitName]), analyzeMuTauPairs)
        optRESnSVfit['sequence'] += analyzeMuTauPairs
    
        analyzeDiTauPairs = cms.EDAnalyzer("ToyMCnSVfitDiTauPairResolutionAnalyzer",
            srcNSVfitEventHypothesis = cms.InputTag(composeModuleName(["nSVfitDiTauPairHypotheses", optRESnSVfitName, optLLnSVfitName])),
            srcGenParticles = cms.InputTag('genTausFromZs'),
            dqmDirectory = cms.string("/%s/nSVfitDiTauPairs/%s" % (optRESnSVfitName, optLLnSVfitName))
        )
        setattr(process, composeModuleName(["analyzeDiTauPairs", optRESnSVfitName]), analyzeDiTauPairs)
        optRESnSVfit['sequence'] += analyzeDiTauPairs
    
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('plotsNSVfitToyMCresolution.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/data2/veelken/CMSSW_4_1_x/skims/ZtoMuTau/DYtautau_spring11_powhegZ2_1_1_XvY.root'
    )
)

process.printGenParticleList = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint = cms.untracked.int32(100)
)

process.o = cms.Path(
    #process.printGenParticleList
    process.produceToyParticles
)

# before starting to process 1st event, print event content
process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")
process.filterFirstEvent = cms.EDFilter("EventCountFilter",
    numEvents = cms.int32(1)
)
process.p1 = cms.Path(process.filterFirstEvent + process.printEventContent)
process.p2 = cms.Path(process.runIdealResolution)
process.p3 = cms.Path(process.runRealisticResolution)

process.q = cms.EndPath(process.savePlots)

process.s = cms.Schedule(process.o, process.p1, process.p2, process.p3, process.q)

# print-out all python configuration parameter information
#print process.dumpPython()
