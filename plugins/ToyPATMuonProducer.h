#ifndef TauAnalysis_GenSimTools_ToyPATMuonProducer_h
#define TauAnalysis_GenSimTools_ToyPATMuonProducer_h

/** \class ToyPATMuonProducer
 *
 * Produce "toy MC" pat::Muon objects 
 * by smearing Monte Carlo "truth" information 
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: ToyPATElectronProducer.h,v 1.3 2010/09/28 11:23:34 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/GenSimTools/interface/ToyParticleProducerBase.h"

class ToyPATMuonProducer : public ToyParticleProducerBase
{
 public:
  // constructor 
  explicit ToyPATMuonProducer(const edm::ParameterSet&);
  ~ToyPATMuonProducer();
  
  void produce(edm::Event&, const edm::EventSetup&);
};

#endif
