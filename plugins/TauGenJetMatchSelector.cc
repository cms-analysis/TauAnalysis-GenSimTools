#include "TauAnalysis/GenSimTools/plugins/TauGenJetMatchSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <vector>

TauGenJetMatchSelector::TauGenJetMatchSelector(const edm::ParameterSet& cfg)
{ 
  //std::cout << "<TauGenJetMatchSelector::TauGenJetMatchSelector>:" << std::endl;

  srcGenTauLeptons_ = cfg.getParameter<edm::InputTag>("srcGenTauLeptons");
  //std::cout << " src = " << srcGenTauLeptons_ << std::endl;
  
  srcGenParticles_ = cfg.getParameter<edm::InputTag>("srcGenParticles");
  //std::cout << " srcGenParticles = " << srcGenParticles_ << std::endl;

  srcTauGenJets_ = cfg.getParameter<edm::InputTag>("srcTauGenJets");
  //std::cout << " srcTauGenJets = " << srcTauGenJets_ << std::endl;
  
  dRmatchGenParticle_ = cfg.getParameter<double>("dRmatchGenParticle");
  //std::cout << " dRmatchGenParticle = " << dRmatchGenParticle_ << std::endl;
  
  dRmatchTauGenJet_ = cfg.getParameter<double>("dRmatchTauGenJet");
  //std::cout << " dRmatchTauGenJet = " << dRmatchTauGenJet_ << std::endl;
  
  produces<reco::GenJetCollection>("");
}

TauGenJetMatchSelector::~TauGenJetMatchSelector()
{
// nothing to be done yet...
}

void TauGenJetMatchSelector::produce(edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<TauGenJetMatchSelector::produce>:" << std::endl;

  edm::Handle<reco::GenParticleCollection> genTauLeptons;
  evt.getByLabel(srcGenTauLeptons_, genTauLeptons);
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByLabel(srcGenParticles_, genParticles);
  
  edm::Handle<reco::GenJetCollection> genTauJets;
  evt.getByLabel(srcTauGenJets_, genTauJets);

  std::auto_ptr<reco::GenJetCollection> matchedTauGenJets(new reco::GenJetCollection());

  std::vector<int> tauLeptonPdgIds;
  tauLeptonPdgIds.push_back(-15);
  tauLeptonPdgIds.push_back(+15);

  for ( reco::GenParticleCollection::const_iterator genTauLepton = genTauLeptons->begin();
	genTauLepton != genTauLeptons->end(); ++genTauLepton ) {
    const reco::GenParticle* matchedGenTauLepton = 
      findGenParticle(genTauLepton->p4(), *genParticles, dRmatchGenParticle_, -1, &tauLeptonPdgIds);

    if ( !matchedGenTauLepton ) { 
      edm::LogWarning ("TauGenJetMatchSelector::produce") 
	<< " Failed to match genTauLepton: Pt = " << genTauLepton->pt() << ","
	<< " eta = " << genTauLepton->eta() << "," 
	<< " phi = " << genTauLepton->phi() << " !!" << std::endl;
      std::cout << "genParticles:" << std::endl;
      for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin();
	    genParticle != genParticles->end(); ++genParticle ) {
	std::cout << "gen. Particle (pdgId = " << genParticle->pdgId() << "): Pt = " << genParticle->pt() << ","
		  << " eta = " << genParticle->eta() << ", phi = " << genParticle->phi() << std::endl;
      }
      continue;
    }

//--- follow generator level particle history to last tau lepton entry
//    (right before the tau lepton decays),
//    in order to exclude from the list of tau lepton decay products
//    photons radiated-off the tau lepton
//    (in Z --> tau+ tau- --> tau+ tau- gamma --> ... events)
    const reco::GenParticle* matchedGenTauLepton_bak = 0;
    while ( matchedGenTauLepton_bak != matchedGenTauLepton ) {
      matchedGenTauLepton_bak = matchedGenTauLepton;
      unsigned numGenDaughters = matchedGenTauLepton->numberOfDaughters();
      for ( unsigned iGenDaughter = 0; iGenDaughter < numGenDaughters; ++iGenDaughter ) {
	const reco::GenParticle* genDaughter = matchedGenTauLepton->daughterRef(iGenDaughter).get();

	if ( genDaughter->pdgId() == matchedGenTauLepton->pdgId() ) {
	  matchedGenTauLepton = genDaughter;
	  break;
	}
      }
    }

    reco::Candidate::LorentzVector genTauJetMomentum = getVisMomentum(matchedGenTauLepton, &(*genParticles));

    const reco::GenJet* matchedTauGenJet = 0;
    double dRmin = 1.e+3;
    for ( reco::GenJetCollection::const_iterator genTauJet = genTauJets->begin();
	  genTauJet != genTauJets->end(); ++genTauJet ) {
      double dR = deltaR(genTauJet->p4(), genTauJetMomentum);
      if ( dR < dRmatchTauGenJet_ && dR < dRmin ) {
	matchedTauGenJet = &(*genTauJet);
	dRmin = dR;
      }
    }

    if ( !matchedTauGenJet ) {
      edm::LogWarning ("TauGenJetMatchSelector::produce") 
	<< " Failed to match genTauJet: Pt = " << genTauJetMomentum.pt() << ","
	<< " eta = " << genTauJetMomentum.eta() << "," 
	<< " phi = " << genTauJetMomentum.phi() << " !!" << std::endl;
      std::cout << "gen. tau lepton: Pt = " << genTauLepton->pt() << ","
		<< " eta = " << genTauLepton->eta() << "," 
		<< " phi = " << genTauLepton->phi() << " !!" << std::endl;
      std::vector<const reco::GenParticle*> genStableDecayProducts;
      findDaughters(matchedGenTauLepton, genStableDecayProducts, +1);
      std::cout << "gen. stable decay products:" << std::endl;
      for ( std::vector<const reco::GenParticle*>::const_iterator genDecayProduct = genStableDecayProducts.begin();
	    genDecayProduct != genStableDecayProducts.end(); ++genDecayProduct ) {
	std::cout << " gen. Particle (pdgId = " << (*genDecayProduct)->pdgId() << "): Pt = " << (*genDecayProduct)->pt() << ","
		  << " eta = " << (*genDecayProduct)->eta() << ", phi = " << (*genDecayProduct)->phi() << std::endl;
      }
      for ( reco::GenJetCollection::const_iterator genTauJet = genTauJets->begin();
	    genTauJet != genTauJets->end(); ++genTauJet ) {
	std::cout << "gen. tau-jet (decay mode = " << JetMCTagUtils::genTauDecayMode(*genTauJet) << "):"
		  << " Pt = " << genTauJet->pt() << ", eta = " << genTauJet->eta() << ", phi = " << genTauJet->phi() << std::endl;
	std::cout << "constituents:" << std::endl;
	const std::vector<const reco::GenParticle*> genTauJetConstituents = genTauJet->getGenConstituents();
	for ( std::vector<const reco::GenParticle*>:: const_iterator genConstituent = genTauJetConstituents.begin();
	      genConstituent != genTauJetConstituents.end(); ++genConstituent ) {
	  std::cout << " gen. Particle (pdgId = " << (*genConstituent)->pdgId() << "): Pt = " << (*genConstituent)->pt() << ","
		    << " eta = " << (*genConstituent)->eta() << ", phi = " << (*genConstituent)->phi() << std::endl;
	}
      }
      continue;
    }
    
    matchedTauGenJets->push_back(reco::GenJet(*matchedTauGenJet));
  }

  evt.put(matchedTauGenJets);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TauGenJetMatchSelector);
