#include "TauAnalysis/GenSimTools/plugins/ToyPATMuonProducer.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"

ToyPATMuonProducer::ToyPATMuonProducer(const edm::ParameterSet& cfg)
  : ToyParticleProducerBase(cfg)
{ 
  produces<pat::MuonCollection>("");
}

ToyPATMuonProducer::~ToyPATMuonProducer()
{
// nothing to be done yet...
}

void ToyPATMuonProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<ToyPATMuonProducer::produce>:" << std::endl;

  std::auto_ptr<pat::MuonCollection> patMuons(new pat::MuonCollection());
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByLabel(src_, genParticles);

  unsigned numGenParticles = genParticles->size();
  for ( unsigned idx = 0; idx < numGenParticles; ++idx ) {
    reco::GenParticleRef genParticle(genParticles, idx);
    reco::Muon muon;
    reco::Candidate::LorentzVector p4smeared = getP4smeared(genParticle->p4());
    muon.setP4(p4smeared);  
    muon.setCharge(genParticle->charge());
    pat::Muon patMuon(muon);
    patMuon.setGenLepton(genParticle);
    //std::cout << " genParticle: Pt = " << genParticle->pt() << "," 
    //	        << " eta = " << genParticle->eta() << ", phi = " << genParticle->phi() << std::endl;
    //std::cout << "--> patMuon: Pt = " << patMuon.pt() << "," 
    //	        << " eta = " << patMuon.eta() << ", phi = " << patMuon.phi() << std::endl;
    patMuons->push_back(patMuon);
  }

  evt.put(patMuons);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ToyPATMuonProducer);
