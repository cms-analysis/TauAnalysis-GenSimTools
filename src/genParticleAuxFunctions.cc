#include "TauAnalysis/GenSimTools/interface/genParticleAuxFunctions.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <TMath.h>

reco::GenParticleCollection getStableDecayProducts(const reco::GenParticle& genMother)
{
//--- return set of stable particles associated with decay of particle given as function argument

  reco::GenParticleCollection genDecayProducts;

  if ( genMother.numberOfDaughters() == 0 ) { // particle given as function argument is stable
    genDecayProducts.push_back(genMother);
  } else { // particle given as function argument is unstable
    size_t numDaughters = genMother.numberOfDaughters();
    for ( size_t iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
      const reco::GenParticle* genDaughter = dynamic_cast<const reco::GenParticle*>(genMother.daughter(iDaughter));
//--- call function recursively for each "daughter" particle
      reco::GenParticleCollection genDaughterDecayProducts = getStableDecayProducts(*genDaughter);
      genDecayProducts.insert(genDecayProducts.end(), genDaughterDecayProducts.begin(), genDaughterDecayProducts.end());
    }
  }

  return genDecayProducts;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

reco::Particle::LorentzVector getGenMissingTransverseMomentum(const reco::GenParticleCollection& genParticles)
{
  reco::Particle::LorentzVector genMissingTransverseMomentum(0,0,0,0);
  for ( reco::GenParticleCollection::const_iterator genParticle = genParticles.begin();
	genParticle != genParticles.end(); ++genParticle ) {
    const reco::Particle::LorentzVector& genParticleMomentum = genParticle->p4();
    int pdgId = genParticle->pdgId();

    if ( pdgId == -12 || pdgId == +12 ||
         pdgId == -14 || pdgId == +14 ||
         pdgId == -16 || pdgId == +16 ) {
      genMissingTransverseMomentum += genParticleMomentum;
    }
  }

  return genMissingTransverseMomentum;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

int getMatchingGenParticlePdgId(const reco::Particle::LorentzVector& recoMomentum,
				const reco::GenParticleCollection& genParticleCollection,
				const std::vector<int>* skipPdgIds, bool useStatusTwoParticles)
{
//--- select genParticles matching direction of reconstructed particle
//    within cone of size dR = 0.5;
//    require generated transverse momentum to be at least half of reconstructed transverse momentum
  reco::GenParticleCollection matchingGenParticles;
  for ( reco::GenParticleCollection::const_iterator genParticle = genParticleCollection.begin();
	genParticle != genParticleCollection.end(); ++genParticle ) {

//--- skip "documentation line" entries
//    (copied over to reco::GenParticle from HepMC product)
		if ( genParticle->status() == 3 && useStatusTwoParticles) continue;
//--- exclude status = 2 particles (i.e. taus) so that tau daughter is matched
		if ( genParticle->status() > 1 && !useStatusTwoParticles) continue;


//--- skip "invisible" particles (e.g. neutrinos);
//    configurable via list of pdgIds given as function argument
    bool skip = false;
    if ( skipPdgIds ) {
      for ( std::vector<int>::const_iterator skipPdgId = skipPdgIds->begin();
	    skipPdgId != skipPdgIds->end(); ++skipPdgId ) {
	if ( TMath::Abs(*skipPdgId) == TMath::Abs(genParticle->pdgId()) ) skip = true;
      }
    } else {
      if ( TMath::Abs(genParticle->pdgId()) == 12 ||
	   TMath::Abs(genParticle->pdgId()) == 14 ||
	   TMath::Abs(genParticle->pdgId()) == 16 ) skip = true;
    }
    if ( skip ) continue;

    if ( genParticle->pt() > 0.50*recoMomentum.pt() &&
	 reco::deltaR(genParticle->p4(), recoMomentum) < 0.5 ) {
      matchingGenParticles.push_back(*genParticle);
    }
  }

//--- find highest Pt matching genParticle
  double ptMax = -1.;
  int pdgId = -1;
  for ( reco::GenParticleCollection::const_iterator matchingGenParticle = matchingGenParticles.begin();
	matchingGenParticle != matchingGenParticles.end(); ++matchingGenParticle ) {

    if ( matchingGenParticle->pt() > ptMax ) {
      pdgId = matchingGenParticle->pdgId();
      ptMax = matchingGenParticle->pt();
    }
  }

  return pdgId;
}

const reco::GenParticle* findMotherWithPdgId(const reco::GenParticle* genParticle, unsigned absPdgId) 
{
  // If this doesn't have any valid mother, return null
  const reco::GenParticle* retVal = 0;
  
  if ( genParticle ) {
    // Get list of mothers.  Not sure why this would be more than one.
    size_t nGenMothers = genParticle->numberOfMothers();
    for ( size_t iGenMother = 0; iGenMother < nGenMothers; ++iGenMother ) {
      const reco::GenParticle* genMother = dynamic_cast<const reco::GenParticle*>(genParticle->mother(iGenMother));
      if ( genMother ) {
	unsigned genMotherPdgId = TMath::Abs(genMother->pdgId());
	if ( genMotherPdgId == absPdgId ) {
	  retVal = genMother;
        } else {
	  const reco::GenParticle* genGrandMother = findMotherWithPdgId(genMother, absPdgId);
	  if ( genGrandMother != 0 ) retVal = genGrandMother;
	}
      }

      if ( retVal != 0 ) break;
    }
  }

  return retVal;
}

























