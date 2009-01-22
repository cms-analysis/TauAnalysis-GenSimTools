#ifndef HiggsAnalysis_GenSimTools_GenMETFromNeutralsProducer_
#define HiggsAnalysis_GenSimTools_GenMETFromNeutralsProducer_

// system include files
#include <memory>
#include <string>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"


/**\class GenMETFromNeutralsProducer 
\brief builds a Candidate from all stable invisible particles (neutrinos, neutralinos, etc.) in the event.

\author Michal Bluj
\date   november 2008
*/
class GenMETFromNeutralsProducer : public edm::EDProducer {
 public:

  explicit GenMETFromNeutralsProducer(const edm::ParameterSet&);

  ~GenMETFromNeutralsProducer();
  
  virtual void produce(edm::Event&, const edm::EventSetup&);

  virtual void beginJob(const edm::EventSetup & c);

 private:
   
  /// Input GenCandidates
  edm::InputTag inputTagGenParticles_;

  /// if yes, use only neutrinos from direct taus
  bool useNeutrinosFromTaus_;

  /// kinematic requirements
  /// (ignored when useNeutrinosFromTaus_ == true)
  double maxEta_, minPt_;

  /// if yes, BSM particles will be excluded 
  /// (ignored when useNeutrinosFromTaus_ == true)
  bool excludeBSM_;
  
  /// verbose ?
  bool verbose_;

  std::set<int> neutrinosPdgId_;

};

#endif
