#ifndef TauAnalysis_GenSimTools_GenPhaseSpaceEventInfoProducer_h
#define TauAnalysis_GenSimTools_GenPhaseSpaceEventInfoProducer_h

/** \class GenPhaseSpaceEventInfoProducer
 *
 * Produce event level information about about phase-space simulated in Monte Carlo
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: GenPhaseSpaceEventInfoProducer.h,v 1.2 2010/01/25 17:25:16 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class GenPhaseSpaceEventInfoProducer : public edm::EDProducer 
{
 public:
  // constructor 
  explicit GenPhaseSpaceEventInfoProducer(const edm::ParameterSet&);
  ~GenPhaseSpaceEventInfoProducer();
  
  void produce(edm::Event&, const edm::EventSetup&);
 
 private:
  // source collection label
  edm::InputTag srcGenEventInfo_;
  edm::InputTag srcGenParticles_;
};

#endif
