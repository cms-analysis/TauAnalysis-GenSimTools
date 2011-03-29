#ifndef TauAnalysis_GenSimTools_ToyPATMEtProducer_h
#define TauAnalysis_GenSimTools_ToyPATMEtProducer_h

/** \class ToyPATMEtProducer
 *
 * Produce "toy MC" pat::MET objects 
 * by smearing Monte Carlo "truth" information 
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

#include "TauAnalysis/GenSimTools/interface/ToyParticleProducerBase.h"

#include <TRandom3.h>

class ToyPATMEtProducer : public edm::EDProducer
{
 public:
  explicit ToyPATMEtProducer(const edm::ParameterSet&);
  ~ToyPATMEtProducer();
  
  void produce(edm::Event&, const edm::EventSetup&);

 private:
  edm::InputTag srcGenParticles_;

  edm::InputTag srcElectrons_;
  edm::InputTag srcMuons_;
  edm::InputTag srcTaus_;

  double resolutionX_;
  double resolutionY_;

  edm::InputTag srcGenMEt_;

  TRandom3 rnd_;
};

#endif
