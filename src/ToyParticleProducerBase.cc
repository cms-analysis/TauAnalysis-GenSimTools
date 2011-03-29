#include "TauAnalysis/GenSimTools/interface/ToyParticleProducerBase.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

ToyParticleProducerBase::ToyParticleProducerBase(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src");

  resolution_ = cfg.getParameter<double>("resolution");
}

ToyParticleProducerBase::~ToyParticleProducerBase()
{
// nothing to be done yet...
}

reco::Candidate::LorentzVector ToyParticleProducerBase::getP4smeared(const reco::Candidate::LorentzVector& p4gen)
{
  double rndSmearFactor = rnd_.Gaus(1., resolution_);
  reco::Candidate::LorentzVector p4smeared(p4gen);
  p4smeared *= rndSmearFactor;
  return p4smeared;
}
 
