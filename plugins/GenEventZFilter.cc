

// system include files
#include <memory>


#include "TauAnalysis/GenSimTools/plugins/GenEventZFilter.h"


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TauAnalysis/CandidateTools/interface/FetchCollection.h"
#include "DataFormats/Common/interface/RefToPtr.h"


#include "PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace edm;
using namespace reco;
using namespace GenParticlesHelper;

GenEventZFilter::GenEventZFilter(const edm::ParameterSet& iConfig)
  : 
  inputTagGenParticles_ ( iConfig.getParameter<InputTag>("GenParticles") ),
  verbose_ ( iConfig.getUntrackedParameter<bool>("verbose", 
						 true) )
{ 
  
  vector<int> dIDs = iConfig.getParameter< vector<int> >("Z0DaughtersPdgId");
  Z0DaughtersPdgId_.insert( dIDs.begin(), dIDs.end() );
  
}


GenEventZFilter::~GenEventZFilter() {}



bool
GenEventZFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;


  Handle<GenParticleCollection> genParticles;
  pf::fetchCollection(genParticles, 
		      inputTagGenParticles_, 
		      iEvent );


  if(verbose_) cout<<"GenEventZFilter"<<endl;

  GenParticleRefVector allStatus3Z0s;  
  findParticles( *genParticles,
		 allStatus3Z0s, 23, 3);
  

  for( IGR iZ0=allStatus3Z0s.begin(); iZ0!=allStatus3Z0s.end(); ++iZ0) {
    if(verbose_) cout<<(*iZ0)<<endl;

    const GenParticleRefVector& daughterRefs = (*iZ0)->daughterRefVector();
    
    assert( !daughterRefs.empty() );
    
    int daughterPdgId = abs( (*daughterRefs.begin())->pdgId() );

    if(verbose_) cout<<"\tdaughter pdgID = "<<daughterPdgId<<endl;
        
    if(  Z0DaughtersPdgId_.find( daughterPdgId ) != Z0DaughtersPdgId_.end() ) {
      if(verbose_) cout<<"\tPASSED"<<endl;
      return true;
    }
  }

  return false;

}


//define this as a plug-in
DEFINE_FWK_MODULE(GenEventZFilter);
