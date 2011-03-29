#ifndef TauAnalysis_GenSimTools_ToyMCnSVfitT1T2ResolutionAnalyzer_h
#define TauAnalysis_GenSimTools_ToyMCnSVfitT1T2ResolutionAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMServices/Core/interface/MonitorElement.h"

#include <string>

template<typename T1, typename T2>
class ToyMCnSVfitT1T2ResolutionAnalyzer : public edm::EDAnalyzer 
{
 public:
  explicit ToyMCnSVfitT1T2ResolutionAnalyzer(const edm::ParameterSet&);
  virtual ~ToyMCnSVfitT1T2ResolutionAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);

  edm::InputTag srcNSVfitEventHypothesis_;
  edm::InputTag srcGenParticles_;

  std::string dqmDirectory_;
  
  MonitorElement* hGenLeg1Pt_;
  MonitorElement* hRecLeg1Pt_;
  MonitorElement* hDeltaLeg1Pt_;

  MonitorElement* hGenLeg2Pt_;
  MonitorElement* hRecLeg2Pt_;
  MonitorElement* hDeltaLeg2Pt_;
  
  MonitorElement* hDPhi_;
  
  MonitorElement* hGenMEt_;
  MonitorElement* hRecMEt_;
  MonitorElement* hDeltaMEt_;

  MonitorElement* hIsValidSolution_;

  MonitorElement* hGenMass_;
  MonitorElement* hRecMassMedian_;
  MonitorElement* hDeltaMassMedian_;
  MonitorElement* hDeltaMassMedianS_;
  MonitorElement* hRecMassMax_;
  MonitorElement* hDeltaMassMax_;
  MonitorElement* hDeltaMassMaxS_;
};

#endif   

