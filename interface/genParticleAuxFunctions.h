#ifndef TauAnalysis_GenSimTools_genParticleAuxFunctions_h
#define TauAnalysis_GenSimTools_genParticleAuxFunctions_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"

reco::GenParticleCollection getStableDecayProducts(const reco::GenParticle& genMother);

reco::Particle::LorentzVector getGenMissingTransverseMomentum(const reco::GenParticleCollection& genParticles);

int getMatchingGenParticlePdgId(const reco::Particle::LorentzVector&, const reco::GenParticleCollection&,
				const std::vector<int>* = 0, bool useStatusTwoParticles = false);

const reco::GenParticle* findMotherWithPdgId(const reco::GenParticle* input, unsigned absPdgId);

#endif
