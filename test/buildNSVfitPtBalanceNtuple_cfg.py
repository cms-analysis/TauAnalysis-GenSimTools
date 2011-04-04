import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(int(1e6)))

process.load("Configuration.GenProduction.PYTHIA6_Tauola_gg_H100_tautau_7TeV_cff")

prodMech = 'bb'
if prodMech == 'bb':
    toModify = process.generator.PythiaParameters.processParameters
    assert(toModify[1].startswith('MSUB(157)'))
    toModify[1] = 'MSUB(186)= 1   ! gg->QQbarH (MSSM)'
    toModify.append(
        'KFPR(186,2)= 5 ! Q = b Registered by Alexandre.Nikitenko@cern.ch',
    )

process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("Configuration.StandardSequences.Services_cff")

process.load("TauAnalysis.GenSimTools.gen_decaysFromAHs_cfi")
process.load("RecoTauTag.TauTagTools.TauTruthProduction_cfi")

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("ptBalance_ntuple.root")
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
    resonanceMass = cms.double(100),
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
