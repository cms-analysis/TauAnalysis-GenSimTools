

// system include files
#include <memory>

#include "TauAnalysis/GenSimTools/plugins/GenTauPairDecayFilter.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "TauAnalysis/GenSimTools/plugins/GenTauJetsHelper.h"

#include "TauAnalysis/CandidateTools/interface/FetchCollection.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace edm;
using namespace reco;

GenTauPairDecayFilter::GenTauPairDecayFilter(const edm::ParameterSet& iConfig)
  : 
  inputTagTauGenJets_ ( iConfig.getParameter<InputTag>("tauGenJets") ),
  verbose_ ( iConfig.getUntrackedParameter<bool>("verbose", 
						 true) )
{ 
  
  vector<int> dID1s 
    = iConfig.getParameter< vector<int> >("tau1DaughtersPdgId");
  tau1DaughtersPdgId_.insert( dID1s.begin(), dID1s.end() );
  vector<int> dID2s 
    = iConfig.getParameter< vector<int> >("tau2DaughtersPdgId");
  tau2DaughtersPdgId_.insert( dID2s.begin(), dID2s.end() );
  
}


GenTauPairDecayFilter::~GenTauPairDecayFilter() {}



bool
GenTauPairDecayFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;


  Handle<GenJetCollection> genJets;
  pf::fetchCollection(genJets, 
		      inputTagTauGenJets_, 
		      iEvent );


  if(verbose_) cout<<"GenTauPairDecayFilter"<<endl;
  
  // need 2 taus 
  if( genJets->size()!=2 ) return false;

  reco::GenJetCollection::const_iterator iTau=genJets->begin();
  const GenJet& tau1 = *iTau; ++iTau;
  const GenJet& tau2 = *iTau;

  bool result1 = GenTauJetsHelper::goodTau( tau1, tau1DaughtersPdgId_) && GenTauJetsHelper::goodTau( tau2, tau2DaughtersPdgId_ );  
  bool result2 = GenTauJetsHelper::goodTau( tau1, tau2DaughtersPdgId_) && GenTauJetsHelper::goodTau( tau2, tau1DaughtersPdgId_ );  

  if(verbose_) {
    std::cout<<"Tau1 type "<<tau1.print()<<std::endl
	     <<"Tau2 type "<<tau2.print()<<std::endl; 
    std::cout<<"  result: "<<result1<<" OR "<<result2<<"="<<(result1||result2)<<std::endl;
  }
       
  return ( result1 || result2 );
}


//define this as a plug-in
DEFINE_FWK_MODULE(GenTauPairDecayFilter);
