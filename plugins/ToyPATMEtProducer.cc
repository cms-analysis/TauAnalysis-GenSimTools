#include "TauAnalysis/GenSimTools/plugins/ToyPATMEtProducer.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <TMath.h>

ToyPATMEtProducer::ToyPATMEtProducer(const edm::ParameterSet& cfg)
{ 
  if ( cfg.exists("srcGenMEt") ) {
    srcGenMEt_ = cfg.getParameter<edm::InputTag>("srcGenMEt");
  } else {
    srcGenParticles_ = cfg.getParameter<edm::InputTag>("srcGenParticles");
  }

  if ( cfg.exists("srcElectrons") ) srcElectrons_ = cfg.getParameter<edm::InputTag>("srcElectrons");
  if ( cfg.exists("srcMuons")     ) srcMuons_     = cfg.getParameter<edm::InputTag>("srcMuons"); 
  if ( cfg.exists("srcTaus")      ) srcTaus_      = cfg.getParameter<edm::InputTag>("srcTaus");

  resolutionX_ = cfg.getParameter<double>("resolutionX");
  resolutionY_ = cfg.getParameter<double>("resolutionY");

  produces<pat::METCollection>("");
}

ToyPATMEtProducer::~ToyPATMEtProducer()
{
// nothing to be done yet...
}

void ToyPATMEtProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<ToyPATMEtProducer::produce>:" << std::endl;

  reco::Candidate::LorentzVector genMEtP4;
  if ( srcGenMEt_.label() != "" ) {
    typedef edm::View<reco::MET> METView;
    edm::Handle<METView> genMEt;
    evt.getByLabel(srcGenMEt_, genMEt);
    if ( !genMEt->size() == 1 ) 
      throw cms::Exception("ToyPATMEtProducer::produce") 
	<< "Failed to find unique gen. MET object !!\n";
    genMEtP4 = genMEt->front().p4();
  } else {
    edm::Handle<reco::GenParticleCollection> genParticles;
    evt.getByLabel(srcGenParticles_, genParticles);
    for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin();
	  genParticle != genParticles->end(); ++genParticle ) {
      if ( genParticle->status() == 1 && isNeutrino(&(*genParticle)) ) genMEtP4 += genParticle->p4();
    }
  }
  
  double metPx = rnd_.Gaus(0., resolutionX_);
  double metPy = rnd_.Gaus(0., resolutionY_);

  if ( srcElectrons_.label() != "" ) {
    edm::Handle<pat::ElectronCollection> patElectrons;
    evt.getByLabel(srcElectrons_, patElectrons);
    for ( pat::ElectronCollection::const_iterator patElectron = patElectrons->begin();
	  patElectron != patElectrons->end(); ++patElectron ) {
      if ( patElectron->genLepton() ) {
	reco::Candidate::LorentzVector genElectronMomentum = patElectron->genLepton()->p4();
	metPx -= (patElectron->px() - genElectronMomentum.px());
	metPy -= (patElectron->py() - genElectronMomentum.py());
      }
    }
  }

  if ( srcMuons_.label() != "" ) {
    edm::Handle<pat::MuonCollection> patMuons;
    evt.getByLabel(srcMuons_, patMuons);
    for ( pat::MuonCollection::const_iterator patMuon = patMuons->begin();
	  patMuon != patMuons->end(); ++patMuon ) {
      if ( patMuon->genLepton() ) { 
	reco::Candidate::LorentzVector genMuonMomentum = patMuon->genLepton()->p4();
	metPx -= (patMuon->px() - genMuonMomentum.px());
	metPy -= (patMuon->py() - genMuonMomentum.py());
      }
    }
  }
  
  if ( srcTaus_.label() != "" ) {
    edm::Handle<pat::TauCollection> patTaus;
    evt.getByLabel(srcTaus_, patTaus);
    for ( pat::TauCollection::const_iterator patTau = patTaus->begin();
	  patTau != patTaus->end(); ++patTau ) {
      if ( patTau->genJet() ) {
	reco::Candidate::LorentzVector genTauMomentum = getVisMomentum(patTau->genJet()->getGenConstituents());
	metPx -= (patTau->px() - genTauMomentum.px());
	metPy -= (patTau->py() - genTauMomentum.py());
      }
    }
  }
  
  std::auto_ptr<pat::METCollection> patMETs(new pat::METCollection());
  reco::MET met;
  met.setP4(reco::Candidate::LorentzVector(metPx, metPy, 0., TMath::Sqrt(metPx*metPx + metPy*metPy)));
  reco::GenMET genMEt;
  genMEt.setP4(genMEtP4);
  pat::MET patMEt(met);
  patMEt.setGenMET(genMEt);
  patMETs->push_back(patMEt);

  evt.put(patMETs);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ToyPATMEtProducer);
