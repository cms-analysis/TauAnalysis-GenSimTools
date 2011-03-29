#include "TauAnalysis/GenSimTools/plugins/ToyPATElectronProducer.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

ToyPATElectronProducer::ToyPATElectronProducer(const edm::ParameterSet& cfg)
  : ToyParticleProducerBase(cfg)
{ 
  produces<pat::ElectronCollection>("");
}

ToyPATElectronProducer::~ToyPATElectronProducer()
{
// nothing to be done yet...
}

void ToyPATElectronProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<ToyPATElectronProducer::produce>:" << std::endl;

  std::auto_ptr<pat::ElectronCollection> patElectrons(new pat::ElectronCollection());

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByLabel(src_, genParticles);

  unsigned numGenParticles = genParticles->size();
  for ( unsigned idx = 0; idx < numGenParticles; ++idx ) {
    reco::GenParticleRef genParticle(genParticles, idx);
    reco::GsfElectron gsfElectron;
    reco::Candidate::LorentzVector p4smeared = getP4smeared(genParticle->p4());
    gsfElectron.setP4(p4smeared);  
    gsfElectron.setCharge(genParticle->charge());
    pat::Electron patElectron(gsfElectron);
    patElectron.setGenLepton(genParticle);
    //std::cout << " genParticle: Pt = " << genParticle->pt() << "," 
    //	        << " eta = " << genParticle->eta() << ", phi = " << genParticle->phi() << std::endl;
    //std::cout << "--> patElectron: Pt = " << patElectron.pt() << "," 
    //	        << " eta = " << patElectron.eta() << ", phi = " << patElectron.phi() << std::endl;
    patElectrons->push_back(patElectron);
  }

  evt.put(patElectrons);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ToyPATElectronProducer);
