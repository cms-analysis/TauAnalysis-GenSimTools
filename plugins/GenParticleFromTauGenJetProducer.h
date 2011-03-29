#ifndef TauAnalysis_GenSimTools_GenParticleFromTauGenJetProducer_h
#define TauAnalysis_GenSimTools_GenParticleFromTauGenJetProducer_h

/** \class ToyPATMEtProducer
 *
 * Auxiliary class to "convert" constituents of TauGenJet
 * in GenParticleCollection format
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: ToyPATMEtProducer.h,v 1.3 2010/09/28 11:23:34 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class GenParticleFromTauGenJetProducer : public edm::EDProducer
{
 public:
  // constructor 
  explicit GenParticleFromTauGenJetProducer(const edm::ParameterSet&);
  ~GenParticleFromTauGenJetProducer();
  
  void produce(edm::Event&, const edm::EventSetup&);

 private:
  edm::InputTag src_;

  StringCutObjectSelector<reco::GenParticle>* cut_;
};

#endif
