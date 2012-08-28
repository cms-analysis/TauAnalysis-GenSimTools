#ifndef TauAnalysis_GenSimTools_GenVertexProducer_h
#define TauAnalysis_GenSimTools_GenVertexProducer_h

/** \class GenVertexProducer
 *
 * Auxiliary class to find Monte Carlo truth event vertex of hard-scatter interaction.
 * The generated event vertex is taken as the production vertex of the first genParticle
 * matching any of the pdgIds specified via configuration parameters.
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: GenVertexProducer.h,v 1.4 2012/04/07 15:44:43 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>

class GenVertexProducer : public edm::EDProducer 
{
 public:
  explicit GenVertexProducer(const edm::ParameterSet&);
  ~GenVertexProducer();
  
  void produce(edm::Event&, const edm::EventSetup&);
  
 private:
  edm::InputTag srcGenParticles_;

  typedef std::vector<int> vint;
  vint pdgIds_;
};

#endif
