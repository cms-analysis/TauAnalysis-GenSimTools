#ifndef __HiggsAnalysis_GenSimTools_GenEventZFilter__
#define __HiggsAnalysis_GenSimTools_GenEventZFilter__

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include <stdio.h>
#include <set>

class GenEventZFilter : public edm::EDFilter {
 public:
  typedef reco::GenParticleRefVector::const_iterator IGR;

  explicit GenEventZFilter(const edm::ParameterSet&);
  ~GenEventZFilter();


 private:
  virtual void beginJob(const edm::EventSetup&) {}
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() {}
    

  // ----------member data ---------------------------

  edm::InputTag     inputTagGenParticles_;
  std::set<int>     Z0DaughtersPdgId_;

  bool              verbose_;


};

#endif

