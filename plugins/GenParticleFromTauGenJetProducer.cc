#include "TauAnalysis/GenSimTools/plugins/GenParticleFromTauGenJetProducer.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include <string>

GenParticleFromTauGenJetProducer::GenParticleFromTauGenJetProducer(const edm::ParameterSet& cfg)
  : cut_(0)
{ 
  src_ = cfg.getParameter<edm::InputTag>("src");

  if ( cfg.exists("cut") ) {
    std::string cut_string = cfg.getParameter<std::string>("cut");
    cut_ = new StringCutObjectSelector<reco::GenParticle>(cut_string);
  }

  produces<reco::GenParticleCollection>("");
}

GenParticleFromTauGenJetProducer::~GenParticleFromTauGenJetProducer()
{
  delete cut_;
}

void GenParticleFromTauGenJetProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<GenParticleFromTauGenJetProducer::produce>:" << std::endl;

  std::auto_ptr<reco::GenParticleCollection> genParticles(new reco::GenParticleCollection());

  edm::Handle<reco::GenJetCollection> genTauJets;
  evt.getByLabel(src_, genTauJets);

  for ( reco::GenJetCollection::const_iterator genTauJet = genTauJets->begin();
	genTauJet != genTauJets->end(); ++genTauJet ) {
    std::vector<const reco::GenParticle*> genDecayProducts = genTauJet->getGenConstituents();
    for ( std::vector<const reco::GenParticle*>::const_iterator genDecayProduct = genDecayProducts.begin();
	  genDecayProduct != genDecayProducts.end(); ++genDecayProduct ) {
      if ( (!cut_) || (*cut_)(**genDecayProduct) ) genParticles->push_back(**genDecayProduct);
    }
  }

  evt.put(genParticles);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GenParticleFromTauGenJetProducer);
