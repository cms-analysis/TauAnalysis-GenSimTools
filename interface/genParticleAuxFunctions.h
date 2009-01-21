#ifndef TauAnalysis_GenSimTools_genParticleAuxFunctions_h
#define TauAnalysis_GenSimTools_genParticleAuxFunctions_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"

reco::GenParticleCollection getStableDecayProducts(edm::Handle<edm::View<reco::GenParticle> >& genParticles, 
						   const reco::GenParticle& genMother);

reco::Particle::LorentzVector getGenVisibleTauMomentum(const reco::CompositePtrCandidate& genTauJet);

reco::Particle::LorentzVector getGenMissingTransverseMomentum(edm::Handle<edm::View<reco::GenParticle> >& genParticles);

#endif
