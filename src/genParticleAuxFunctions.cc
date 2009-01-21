#include "TauAnalysis/GenSimTools/interface/genParticleAuxFunctions.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/TriggerNames.h"

reco::GenParticleCollection getStableDecayProducts(edm::Handle<edm::View<reco::GenParticle> >& genParticles, 
						   const reco::GenParticle& genMother)
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
      reco::GenParticleCollection genDaughterDecayProducts = getStableDecayProducts(genParticles, *genDaughter);
      genDecayProducts.insert(genDecayProducts.end(), genDaughterDecayProducts.begin(), genDaughterDecayProducts.end());
    }
  }

  return genDecayProducts;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

reco::Particle::LorentzVector getGenVisibleTauMomentum(const reco::CompositePtrCandidate& genTauJet)
{
  reco::Particle::LorentzVector genVisibleTauMomentum(0,0,0,0);

  const reco::CompositePtrCandidate::daughters& genDecayProducts = genTauJet.daughterPtrVector();
  for ( reco::CompositePtrCandidate::daughters::const_iterator genDecayProduct = genDecayProducts.begin();
 	genDecayProduct != genDecayProducts.end(); ++genDecayProduct ) {

//--- exclude electron, muon and tau neutrinos from visible energy/momentum calculation
    if ( !((*genDecayProduct)->pdgId() == -12 || (*genDecayProduct)->pdgId() == +12 ||
	   (*genDecayProduct)->pdgId() == -14 || (*genDecayProduct)->pdgId() == +14 ||
	   (*genDecayProduct)->pdgId() == -16 || (*genDecayProduct)->pdgId() == +16) ) genVisibleTauMomentum += (*genDecayProduct)->p4();
  }

  return genVisibleTauMomentum;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

reco::Particle::LorentzVector getGenMissingTransverseMomentum(edm::Handle<edm::View<reco::GenParticle> >& genParticles)
{
  reco::Particle::LorentzVector genMissingTransverseMomentum(0,0,0,0);
  for ( edm::View<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); 
	genParticle != genParticles->end(); ++genParticle ) {
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
