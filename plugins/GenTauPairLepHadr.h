#ifndef __HiggsAnalysis_GenSimTools_GenTauPairLepHadr__
#define __HiggsAnalysis_GenSimTools_GenTauPairLepHadr__

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include <stdio.h>
#include <set>

class GenTauPairLepHadr : public edm::EDFilter {
 public:

  explicit GenTauPairLepHadr(const edm::ParameterSet&);
  ~GenTauPairLepHadr();


 private:
  virtual void beginJob(const edm::EventSetup&) {}
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() {}
    
  /// returns -1 if the tau is leptonic, 1 otherwise
  int tauType(const reco::GenJet& jet) const;

  // ----------member data ---------------------------

  edm::InputTag     inputTagTauGenJets_;
  std::set<int>     tauLeptonicDaughtersPdgId_;

  bool              verbose_;


};

#endif

