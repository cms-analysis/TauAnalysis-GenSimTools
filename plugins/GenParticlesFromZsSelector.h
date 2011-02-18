#ifndef TauAnalysis_GenSimTools_GenParticlesFromZsSelector_h
#define TauAnalysis_GenSimTools_GenParticlesFromZsSelector_h

/** \class GenParticlesFromZsSelector
 *
 * Select tau leptons, muons and electrons
 * produced in Z/gamma* --> l+ l- processes
 *
 * NOTE: Can handle case that virtual Z/gamma is missing in GenEVT record
 * 
 * \author Christian Veelken;
 *
 * \version $Revision: 1.2 $
 *
 * $Id: GenParticlesFromZsSelector.h,v 1.2 2010/09/28 11:23:35 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>

class GenParticlesFromZsSelector : public edm::EDProducer 
{
 public:

  explicit GenParticlesFromZsSelector(const edm::ParameterSet&);

  ~GenParticlesFromZsSelector();
  
  void produce(edm::Event&, const edm::EventSetup&);

 private:

  edm::InputTag src_;

  typedef std::vector<int> vint;
  vint pdgIdsMothers_;
  vint pdgIdsDaughters_;

  int maxDaughters_;
  int minDaughters_;
};

#endif
