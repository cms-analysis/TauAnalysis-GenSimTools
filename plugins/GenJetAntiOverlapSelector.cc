#include "TauAnalysis/CandidateTools/interface/ParticleAntiOverlapSelector.h"

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/JetReco/interface/GenJet.h"

typedef ObjectSelector<ParticleAntiOverlapSelector<reco::GenJet> > GenJetAntiOverlapSelector;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GenJetAntiOverlapSelector);
