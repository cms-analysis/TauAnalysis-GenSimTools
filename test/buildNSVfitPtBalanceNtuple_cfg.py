import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("ANA")
process.load("Configuration.GenProduction.PYTHIA6_Tauola_gg_H100_tautau_7TeV_cff")
print "Forcing mu-had decays"
process.generator.ExternalDecays.Tauola.InputCards.mdtau = 216

options = VarParsing.VarParsing ('analysis')
options.register(
    'bbMode', 0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,
    "bbMode=1 specifies bb associated production")

options.register(
    'mass', 90,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,
    "Mass of A0")

options.parseArguments()

if options.maxEvents < 0:
    options.maxEvents = int(1e5)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(
    options.maxEvents))

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

pythia_config = process.generator.PythiaParameters.processParameters
# Set Mass
assert(pythia_config[4].startswith('RMSS(19)'))
pythia_config[4] = 'RMSS(19)= %0.0f ! m_A' % options.mass
print "Setting mass of A0 to %0.0f" % options.mass

if options.bbMode:
    assert(pythia_config[1].startswith('MSUB(157)'))
    pythia_config[1] = 'MSUB(186)= 1   ! gg->QQbarH (MSSM)'
    pythia_config.append(
        'KFPR(186,2)= 5 ! Q = b Registered by Alexandre.Nikitenko@cern.ch',
    )

process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("Configuration.StandardSequences.Services_cff")

process.load("TauAnalysis.GenSimTools.gen_decaysFromAHs_cfi")
process.load("RecoTauTag.TauTagTools.TauTruthProduction_cfi")

prod_label = options.bbMode and "bbH" or "ggAH"

output_file_name = "ptBalanceData_mass_%i_%s.root" % (options.mass,
                                                      prod_label)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(output_file_name)
)

# Turn off the cuts, we apply them in the ntuple
process.genMuonsFromAHtautauDecaysWithinAcceptance.cut = ""
process.genElectronsFromAHtautauDecaysWithinAcceptance.cut = ""
process.genHadronsFromAHtautauDecaysWithinAcceptance.cut = ""

# Make a matching between our visible decays and mother taus
process.load("RecoTauTag.TauTagTools.RecoTauTruthMatching_cfi")

process.matchMuonicDecay = process.recoTauTruthMatcher.clone(
    src = cms.InputTag("genTausFromAHs"),
    matched = cms.InputTag("genMuonsFromAHtautauDecaysWithinAcceptance"),
)

process.matchHadronicDecay = process.recoTauTruthMatcher.clone(
    src = cms.InputTag("genTausFromAHs"),
    matched = cms.InputTag("genHadronsFromAHtautauDecaysWithinAcceptance"),
)

process.requireAtLeast1genMuonsFromAHtautauDecaysWithinAcceptance = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("genMuonsFromAHtautauDecaysWithinAcceptance"),
    minNumber = cms.uint32(1)
)

process.requireAtLeast1genHadronsFromAHtautauDecaysWithinAcceptance = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("genHadronsFromAHtautauDecaysWithinAcceptance"),
    minNumber = cms.uint32(1)
)

# Find the true tau that is associated to the muonic decay (leg 1)
process.matchTauToVisTauMuon = cms.EDFilter(
    "CandViewGenJetMatchRefSelector",
    src = cms.InputTag("genTausFromAHs"),
    matching = cms.InputTag("matchMuonicDecay"),
    filter = cms.bool(True)
)

# Find the true tau that is associated to the hadronic decay (leg 1)
process.matchTauToVisTauHadron = process.matchTauToVisTauMuon.clone(
    matching = cms.InputTag("matchHadronicDecay")
)

process.makePtBalanceNtuple = cms.EDAnalyzer(
    "NSVfitPtBalanceNtupleProducer",
    leg1Src = cms.InputTag("matchTauToVisTauMuon"),
    leg1VisSrc = cms.InputTag("genMuonsFromAHtautauDecaysWithinAcceptance"),
    leg2Src = cms.InputTag("matchTauToVisTauHadron"),
    leg2VisSrc = cms.InputTag("genHadronsFromAHtautauDecaysWithinAcceptance"),
    resonanceMass = cms.double(options.mass),
)

process.p = cms.Path(
    process.generator
    +process.genParticles
    +process.tauGenJets
    +process.produceGenDecayProductsFromAHs
    +process.matchMuonicDecay
    +process.matchHadronicDecay
    +process.requireAtLeast1genHadronsFromAHtautauDecaysWithinAcceptance
    +process.requireAtLeast1genMuonsFromAHtautauDecaysWithinAcceptance
    +process.matchTauToVisTauMuon
    +process.matchTauToVisTauHadron
    +process.makePtBalanceNtuple
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
