#ifndef TauAnalysis_GenSimTools_ToyParticleProducerBase_h
#define TauAnalysis_GenSimTools_ToyParticleProducerBase_h

/** \class ToyParticleProducerBase
 *
 * Base-class for producing "toy MC" pat::Electron, pat::Muon and pat::Tau objects 
 * by smearing Monte Carlo "truth" information 
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: ToyParticleProducerBase.h,v 1.3 2010/09/28 11:23:34 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Candidate/interface/Candidate.h" 

#include <TRandom3.h>

class ToyParticleProducerBase : public edm::EDProducer 
{
 public:
  // constructor 
  explicit ToyParticleProducerBase(const edm::ParameterSet&);
  ~ToyParticleProducerBase();
  
 protected:
  reco::Candidate::LorentzVector getP4smeared(const reco::Candidate::LorentzVector&);
 
  edm::InputTag src_;

  double resolution_;

  TRandom3 rnd_;
};

#endif
