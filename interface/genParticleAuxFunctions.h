#ifndef TauAnalysis_GenSimTools_genParticleAuxFunctions_h
#define TauAnalysis_GenSimTools_genParticleAuxFunctions_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"

reco::GenParticleCollection getStableDecayProducts(const reco::GenParticle& genMother);

reco::Particle::LorentzVector getGenMissingTransverseMomentum(const reco::GenParticleCollection& genParticles);

reco::Particle::LorentzVector getGenMissingTransverseMomentum(edm::Handle<edm::View<reco::GenParticle> >& genParticles);

int getMatchingGenParticlePdgId(const reco::Particle::LorentzVector&, edm::Handle<reco::GenParticleCollection>&, 
				const std::vector<int>* = 0, bool useStatusTwoParticles = false);

#endif
