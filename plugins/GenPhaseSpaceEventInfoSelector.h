#ifndef TauAnalysis_GenSimTools_GenPhaseSpaceEventInfoSelector_h
#define TauAnalysis_GenSimTools_GenPhaseSpaceEventInfoSelector_h

/** \class GenPhaseSpaceEventInfoSelector
 *
 * Selection of events based on generator level information,
 * in order to avoid overlap in phase-space simulated in different Monte Carlo samples
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: GenPhaseSpaceEventSelector.h,v 1.1 2009/01/23 14:58:12 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "PhysicsTools/Utilities/interface/StringCutObjectSelector.h"

#include "AnalysisDataFormats/TauAnalysis/interface/GenPhaseSpaceEventInfo.h"

#include "PhysicsTools/UtilAlgos/interface/EventSelectorBase.h"

class GenPhaseSpaceEventInfoSelector : public EventSelectorBase
{
 public:
  // constructor 
  explicit GenPhaseSpaceEventInfoSelector(const edm::ParameterSet&);
  ~GenPhaseSpaceEventInfoSelector();
  
  bool operator()(edm::Event&, const edm::EventSetup&);
 
 private:
  // source collection label
  edm::InputTag src_;
 
  StringCutObjectSelector<GenPhaseSpaceEventInfo>* cut_;
};

#endif
