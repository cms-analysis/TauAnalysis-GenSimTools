#include "TauAnalysis/GenSimTools/plugins/GenVertexProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

GenVertexProducer::GenVertexProducer(const edm::ParameterSet& cfg) 
{
  srcGenParticles_ = cfg.getParameter<edm::InputTag>("srcGenParticles");

  pdgIds_ = cfg.getParameter<vint>("pdgIds");

  produces<reco::Vertex>();
}

GenVertexProducer::~GenVertexProducer()
{
// nothing to be done yet...
}

void GenVertexProducer::produce(edm::Event& evt, const edm::EventSetup& es) 
{
  reco::Vertex::Point genVertexPos;
  bool isFound = false;
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByLabel(srcGenParticles_, genParticles);

  for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin();
	genParticle != genParticles->end() && !isFound; ++genParticle ) {
    int genParticle_pdgId = genParticle->pdgId();
    for ( vint::const_iterator pdgId = pdgIds_.begin();
	  pdgId != pdgIds_.end() && !isFound; ++pdgId ) {
      if ( genParticle_pdgId == (*pdgId) ) {
	genVertexPos = genParticle->vertex();
        isFound = true;	
      }
    }
  }

  if ( !isFound ) {
    edm::LogWarning ("GenVertexProducer::produce")
      << "Failed to find gen. Event vertex !!";
    genVertexPos.SetXYZ(0., 0., 0.);
  }

  //std::cout << "<GenVertexProducer::produce>:" << std::endl;
  //std::cout << " genVertex: x = " << genVertexPos.x() << ", y = " << genVertexPos.y() << ", z = " << genVertexPos.z() << std::endl;
  
  std::auto_ptr<reco::Vertex> genVertex(new reco::Vertex(genVertexPos, reco::Vertex::Error()));
  evt.put(genVertex);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GenVertexProducer);
