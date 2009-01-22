
#include "TauAnalysis/GenSimTools/plugins/GenTauJetsHelper.h"


namespace GenTauJetsHelper {

  bool goodTau(const reco::GenJet& jet, const std::set<int>& tauDaughtersPdgId){
      
    std::vector<const reco::GenParticle*> constits = jet.getGenConstituents();
    
    int pdgId = -1;
   
    if( constits.size() == 1 ) {      
      const reco::GenParticle* gen = *(constits.begin());      
      if( abs(gen->pdgId()) == 11 ||
	  abs(gen->pdgId()) == 13 ) {
	
	pdgId = abs(gen->pdgId());
      }
    } 
    else if( constits.size() == 2 ) {
      const reco::GenParticle* gen1 = *(constits.begin());
      const reco::GenParticle* gen2 = *(constits.begin()+1);
      if( abs(gen2->pdgId()) == 22 ) { 
        if( abs(gen1->pdgId()) == 11 ||
            abs(gen1->pdgId()) == 13 )	
	  pdgId = abs(gen1->pdgId());
      }
      else if( abs(gen1->pdgId()) == 22 ) { 
        if( abs(gen2->pdgId()) == 11 ||
            abs(gen2->pdgId()) == 13 )
	  pdgId = abs(gen2->pdgId());
      }

    }
    
    if( tauDaughtersPdgId.find( pdgId ) !=  tauDaughtersPdgId.end() )
      return true;
    else 
      return false;
    
  }

}
