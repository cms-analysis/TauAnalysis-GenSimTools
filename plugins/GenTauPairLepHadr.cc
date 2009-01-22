

// system include files
#include <memory>


#include "TauAnalysis/GenSimTools/plugins/GenTauPairLepHadr.h"


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "TauAnalysis/CandidateTools/interface/FetchCollection.h"



#include <iostream>
#include <vector>

using namespace std;
using namespace edm;
using namespace reco;

GenTauPairLepHadr::GenTauPairLepHadr(const edm::ParameterSet& iConfig)
  : 
  inputTagTauGenJets_ ( iConfig.getParameter<InputTag>("tauGenJets") ),
  verbose_ ( iConfig.getUntrackedParameter<bool>("verbose", 
						 true) )
{ 
  
  vector<int> dIDs 
    = iConfig.getParameter< vector<int> >("tauLeptonicDaughtersPdgId");
  tauLeptonicDaughtersPdgId_.insert( dIDs.begin(), dIDs.end() );
  
}


GenTauPairLepHadr::~GenTauPairLepHadr() {}



bool
GenTauPairLepHadr::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;


  Handle<GenJetCollection> genJets;
  pf::fetchCollection(genJets, 
		      inputTagTauGenJets_, 
		      iEvent );


  if(verbose_) cout<<"GenTauPairLepHadr"<<endl;

  // need 2 taus 
  if( genJets->size()!=2 ) return false;

  typedef reco::GenJetCollection::const_iterator IGJ;

  int result=1;
  for( IGJ iTau=genJets->begin(); iTau!=genJets->end(); ++iTau) {
    const GenJet& tau = *iTau;
    
    int thisTauResult = tauType( tau );  
    result *= thisTauResult;
    if( verbose_ ) 
      cout<<"  result: "<<thisTauResult<<endl;
       
        
  }

  // result > 1 means -1 * -1 or 1 * 1, that is 2 leptonic taus 
  // or 2 hadronic taus. We don't want this. 

  return result<0 ? true : false;
}


int GenTauPairLepHadr::tauType(const reco::GenJet& jet) const {


  vector<const GenParticle*> constits = jet.getGenConstituents();
  typedef vector<const GenParticle*>::const_iterator IC;

  if(verbose_) 
    cout<<"Tau type "<<jet.print()<<endl;

  // cannot be a lepton
  // (neutrinos from tau decay not included in the tau gen jet)
  if( constits.size()!=1 ) return -1;

  // not a lepton
  const GenParticle* gen = *(constits.begin());
  int pdgId = abs(gen->pdgId());
  
  if( pdgId!=11 &&
      pdgId!=13 ) return -1;

  // not the right lepton type
  if( tauLeptonicDaughtersPdgId_.find(pdgId) ==  tauLeptonicDaughtersPdgId_.end() )
    return 0;
    
  // lepton, return pdgId


  return pdgId;
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenTauPairLepHadr);
