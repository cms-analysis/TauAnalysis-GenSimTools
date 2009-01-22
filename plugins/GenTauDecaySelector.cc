#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/JetReco/interface/GenJet.h"

#include "TauAnalysis/GenSimTools/plugins/GenTauJetsHelper.h"

#include <vector>

struct GenTauDecaySelectorDef {

  typedef reco::GenJetCollection collection;
  typedef edm::Handle<collection> HandleToCollection;
  typedef std::vector<const reco::GenJet*> container;
  typedef container::const_iterator const_iterator;

  GenTauDecaySelectorDef(const edm::ParameterSet& iConfig) : 
    verbose_ ( iConfig.getUntrackedParameter<bool>("verbose", true) )
  {
    std::vector<int> dIDs 
      = iConfig.getParameter< std::vector<int> >("tauDaughtersPdgId");
    tauDaughtersPdgId_.insert( dIDs.begin(), dIDs.end() );

  }

  ~GenTauDecaySelectorDef(){}

  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }

  void select( const HandleToCollection & hc, 
               const edm::Event & e,
               const edm::EventSetup& s) {
    
    selected_.clear();
    assert( hc.isValid() ); // collection of GenJets

    for( collection::const_iterator iT = hc->begin(); iT != hc->end(); ++iT ) {
      const reco::GenJet& genJet = *iT;
      if(verbose_) 
	std::cout<<"Tau type "<<genJet.print()<<std::endl;    

      if( GenTauJetsHelper::goodTau(genJet, tauDaughtersPdgId_) )
	selected_.push_back( new reco::GenJet(*iT) );
    }
  }

  size_t size() const { return selected_.size(); }

 private:

  // ----------member data ---------------------------

  container selected_;
  
  std::set<int>     tauDaughtersPdgId_;
  
  bool              verbose_;
  
  
};


typedef ObjectSelector<GenTauDecaySelectorDef> GenTauDecaySelector;

DEFINE_FWK_MODULE( GenTauDecaySelector );
