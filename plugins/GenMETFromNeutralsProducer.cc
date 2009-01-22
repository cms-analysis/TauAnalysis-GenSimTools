#include "TauAnalysis/GenSimTools/plugins/GenMETFromNeutralsProducer.h"


#include "TauAnalysis/CandidateTools/interface/FetchCollection.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"

#include "PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h"

using namespace std;
using namespace edm;
using namespace reco;

GenMETFromNeutralsProducer::GenMETFromNeutralsProducer(const edm::ParameterSet& iConfig) {
  

  inputTagGenParticles_ 
    = iConfig.getParameter<InputTag>("GenParticles");

  useNeutrinosFromTaus_
    = iConfig.getUntrackedParameter<bool>("neutrinosFromTaus",false);

  maxEta_
    = iConfig.getParameter<double>("maxEta");

  minPt_
    = iConfig.getParameter<double>("minPt");

  excludeBSM_ 
    = iConfig.getUntrackedParameter<bool>("excludeBSM",false);

  verbose_ = 
    iConfig.getUntrackedParameter<bool>("verbose",false);

  int smNeutrinosPdgIds[] = { 12, 14, 16 }; 
  int bsmNeutrinosPdgIds[] = { 18, 39,
			       1000022, 2000012, 2000014,
			       2000016, 1000039, 5000039,
			       4000012, 9900012, 9900014,
			       9900016 };
  neutrinosPdgId_.insert(smNeutrinosPdgIds,smNeutrinosPdgIds+3);
  if( !excludeBSM_ && !useNeutrinosFromTaus_ )
    neutrinosPdgId_.insert(bsmNeutrinosPdgIds,bsmNeutrinosPdgIds+12);

  if(verbose_){
    cout<<"GenMETFromNeutralsProducer "<<endl;
    cout<<"\t"<<inputTagGenParticles_<<endl;
    if(useNeutrinosFromTaus_)
      cout<<"\t Only neutrinos from direct taus. "<<endl;
    else {
      cout<<"\t minPt: "<<minPt_<<endl;
      cout<<"\t maxEta: "<<maxEta_<<endl;
      if(excludeBSM_)
	cout<<"\t BSM stable neutral particles excluded. "<<endl;
      else
	cout<<"\t BSM stable neutral particles included. "<<endl;
    }
  }

  produces<METCollection>("");
}



GenMETFromNeutralsProducer::~GenMETFromNeutralsProducer() { }



void GenMETFromNeutralsProducer::beginJob(const edm::EventSetup & es) { }


void GenMETFromNeutralsProducer::produce(Event& iEvent, const EventSetup& iSetup) {
  
  Handle<GenParticleCollection> genParticles;
  pf::fetchCollection(genParticles,
		      inputTagGenParticles_,
                      iEvent );
  const GenParticleCollection& genParts = *genParticles;

  std::auto_ptr<METCollection> 
    pOutSimpleMet(new METCollection() );


  typedef GenParticleCollection::const_iterator IG;
  
  math::XYZTLorentzVector sumInvisMom;
  double sumEt = 0;
  if(useNeutrinosFromTaus_) {
    GenParticleRefVector myTaus;
    GenParticlesHelper::findParticles( genParts, myTaus, 15, 2 );

    for( GenParticlesHelper::IGR iTau = myTaus.begin(); 
	 iTau != myTaus.end(); ++iTau) {

      if( !GenParticlesHelper::isDirect( *iTau ) )
	continue;

      // look for all stable (status 1) daughters
      GenParticleRefVector daughters;
      GenParticlesHelper::findDescendents( *iTau, daughters, 1 );

      for(GenParticlesHelper::IGR iDau = daughters.begin();
	  iDau != daughters.end(); ++iDau ) {
	// kinematic cuts ignored
	if( neutrinosPdgId_.find( abs( (*iDau)->pdgId() ) ) !=  
	    neutrinosPdgId_.end() ) {
	  sumInvisMom += (*iDau)->p4();
	  sumEt += (*iDau)->et();
	}
      }
    }
  }
  else {
    for(IG ig = genParts.begin(); ig != genParts.end(); ++ig ) {
      const GenParticle& gen = *ig;
      if( gen.status() != 1 ) continue; //stable particle

      if( neutrinosPdgId_.find(abs(gen.pdgId())) !=  neutrinosPdgId_.end() &&
	  ( abs(gen.eta()) < maxEta_ || maxEta_ < 0 ) && gen.pt() > minPt_ ) { 
	sumInvisMom += gen.p4();
	sumEt += gen.et();
      }
    }
  }
  
  MET met( sumEt, 
	   math::XYZTLorentzVector( sumInvisMom.px(),
				    sumInvisMom.py(),
				    0,
				    sumInvisMom.pt() ), 
	   math::XYZPoint( 0, 0, 0 )
	   ); 
  if(verbose_) {
    cout<<"GenMETFromNeutralsProducer "<<endl;
    cout<<"\t MET: "<<met.et()<<" GeV"<<endl;
  }

  pOutSimpleMet->push_back( met );
  
  iEvent.put( pOutSimpleMet );
}


DEFINE_FWK_MODULE( GenMETFromNeutralsProducer );
