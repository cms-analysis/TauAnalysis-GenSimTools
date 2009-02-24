#ifndef TauAnalysis_GenSimTools_GenPhaseSpaceEventInfoProducer_h
#define TauAnalysis_GenSimTools_GenPhaseSpaceEventInfoProducer_h

/** \class GenPhaseSpaceEventInfoProducer
 *
 * Produce event level information about about phase-space simulated in Monte Carlo
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: GenPhaseSpaceEventInfoProducer.h,v 1.1 2009/01/23 14:58:12 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

class GenPhaseSpaceEventInfoProducer : public edm::EDProducer 
{
 public:
  // constructor 
  explicit GenPhaseSpaceEventInfoProducer(const edm::ParameterSet&);
  ~GenPhaseSpaceEventInfoProducer();
  
  void produce(edm::Event&, const edm::EventSetup&);
 
 private:
  // source collection label
  edm::InputTag srcGenEventScale_;
  edm::InputTag srcGenParticles_;
};

#endif
