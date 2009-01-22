#ifndef __HiggsAnalysis_GenSimTools_GenTauPairDecayFilter__
#define __HiggsAnalysis_GenSimTools_GenTauPairDecayFilter__

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include <stdio.h>
#include <set>

class GenTauPairDecayFilter : public edm::EDFilter {
 public:

  explicit GenTauPairDecayFilter(const edm::ParameterSet&);
  ~GenTauPairDecayFilter();


 private:
  virtual void beginJob(const edm::EventSetup&) {}
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() {}
    
  // ----------member data ---------------------------

  edm::InputTag     inputTagTauGenJets_;
  std::set<int>     tau1DaughtersPdgId_, tau2DaughtersPdgId_;

  bool              verbose_;


};

#endif

