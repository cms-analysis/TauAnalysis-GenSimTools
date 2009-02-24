#include "TauAnalysis/GenSimTools/plugins/GenPhaseSpaceEventInfoProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "AnalysisDataFormats/TauAnalysis/interface/GenPhaseSpaceEventInfo.h"

GenPhaseSpaceEventInfoProducer::GenPhaseSpaceEventInfoProducer(const edm::ParameterSet& cfg)
{ 
  //std::cout << "<GenPhaseSpaceEventInfoProducer::GenPhaseSpaceEventInfoProducer>:" << std::endl;

  srcGenEventScale_ = cfg.getParameter<edm::InputTag>("srcGenEventScale");
  //std::cout << " srcGenEventScale = " << srcGenEventScale_ << std::endl;

  srcGenParticles_ = cfg.getParameter<edm::InputTag>("srcGenParticles");
  //std::cout << " srcGenParticles = " << srcGenParticles_ << std::endl;

  produces<GenPhaseSpaceEventInfo>("");
}

GenPhaseSpaceEventInfoProducer::~GenPhaseSpaceEventInfoProducer()
{
// nothing to be done yet...
}

void checkGenParticle(const reco::GenParticle& genParticle, int pdgId, reco::Particle::LorentzVector& leadingGenParticleMomentum)
{
//--- skip "documentation line" entries
//    (copied over to reco::GenParticle from HepMC product)
  if ( genParticle.status() == 3 ) return;

//--- require genParticle to be of the specified type
  if ( genParticle.pdgId() == +pdgId || genParticle.pdgId() == -pdgId ) {    

//--- select highest Pt genParticle of matching type
    if ( genParticle.pt() > leadingGenParticleMomentum.pt() ) leadingGenParticleMomentum = genParticle.p4();
  }
}

void GenPhaseSpaceEventInfoProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<GenPhaseSpaceEventInfoProducer::produce>:" << std::endl;

//--- compute Pt(hat)
//
// WARNING: defined only for Monte Carlo samples generated with Pythia
//
  double ptHat = 0.;
  edm::Handle<double> genEventScale;
  evt.getByLabel(srcGenEventScale_, genEventScale);
  if ( genEventScale.isValid() ) {
    ptHat = (*genEventScale);
    //std::cout << "Pt(hat) = " << ptHat << std::endl;
  } else {
    ptHat = -1.;
  }

//--- search generator level particles for highest Pt muon, electron and tau (lepton)
  reco::Particle::LorentzVector leadingGenElectron;
  reco::Particle::LorentzVector leadingGenMuon;
  reco::Particle::LorentzVector leadingGenTauLepton;

  edm::Handle<edm::View<reco::GenParticle> > genParticles;
  evt.getByLabel(srcGenParticles_, genParticles);

  for ( edm::View<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); 
	genParticle != genParticles->end(); ++genParticle ) {

    checkGenParticle(*genParticle, 11, leadingGenElectron);
    checkGenParticle(*genParticle, 13, leadingGenMuon);
    checkGenParticle(*genParticle, 15, leadingGenTauLepton);
  }

//--- create new GenPhaseSpaceEventInfo object
  std::auto_ptr<GenPhaseSpaceEventInfo> genPhaseSpaceEventInfo(new GenPhaseSpaceEventInfo());
  genPhaseSpaceEventInfo->ptHat_ = ptHat;
  genPhaseSpaceEventInfo->leadingGenElectron_ = leadingGenElectron;
  genPhaseSpaceEventInfo->leadingGenMuon_ = leadingGenMuon;
  genPhaseSpaceEventInfo->leadingGenTauLepton_ = leadingGenTauLepton;

//--- add GenPhaseSpaceEventInfo object to the event
  evt.put(genPhaseSpaceEventInfo);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_ANOTHER_FWK_MODULE(GenPhaseSpaceEventInfoProducer);
