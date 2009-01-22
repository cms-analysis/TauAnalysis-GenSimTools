#ifndef HiggsAnalysis_GenSimTools_GenTausJetHelper_h
#define HiggsAnalysis_GenSimTools_GenTausJetHelper_h

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <set>

namespace GenTauJetsHelper {

  typedef std::vector<const reco::GenParticle*>::const_iterator IC;

  /// returns true if the tau decay on requested particle, false otherwise
  bool goodTau(const reco::GenJet& jet, const std::set<int>& tauDaughtersPdgId);

}

#endif
