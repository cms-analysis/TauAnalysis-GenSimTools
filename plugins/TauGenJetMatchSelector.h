#ifndef TauAnalysis_GenSimTools_TauGenJetMatchSelector_h
#define TauAnalysis_GenSimTools_TauGenJetMatchSelector_h

/** \class TauGenJetMatchSelector
 *
 * Produce event level information about about phase-space simulated in Monte Carlo
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: TauGenJetMatchSelector.h,v 1.1 2010/01/20 09:14:08 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class TauGenJetMatchSelector : public edm::EDProducer 
{
 public:
  // constructor 
  explicit TauGenJetMatchSelector(const edm::ParameterSet&);
  ~TauGenJetMatchSelector();
  
  void produce(edm::Event&, const edm::EventSetup&);
 
 private:
  edm::InputTag srcGenTauLeptons_;
  edm::InputTag srcGenParticles_;
  edm::InputTag srcTauGenJets_;
  double dRmatchGenParticle_;
  double dRmatchTauGenJet_;
};

#endif
