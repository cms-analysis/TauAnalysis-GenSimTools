#include "TauAnalysis/GenSimTools/plugins/GenPhaseSpaceEventInfoSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

GenPhaseSpaceEventInfoSelector::GenPhaseSpaceEventInfoSelector(const edm::ParameterSet& cfg)
{ 
  //std::cout << "<GenPhaseSpaceEventInfoSelector::GenPhaseSpaceEventInfoSelector>:" << std::endl;

  src_ = cfg.getParameter<edm::InputTag>("src");
  //std::cout << " src = " << src_ << std::endl;

  std::string cut_string = cfg.getParameter<std::string>("cut");
  if ( cut_string != "" ) {
    cut_ = new StringCutObjectSelector<GenPhaseSpaceEventInfo>(cut_string);
  } else {
    edm::LogWarning ("GenPhaseSpaceEventInfoSelector") << " Empty cut string --> no Events will be rejected !!";
    cut_ = 0;
  }
  //std::cout << " cut_string = " << cut_string << std::endl;
  //std::cout << " cut = " << cut_ << std::endl;
}

GenPhaseSpaceEventInfoSelector::~GenPhaseSpaceEventInfoSelector()
{
  delete cut_;
}

bool GenPhaseSpaceEventInfoSelector::operator()(edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<GenPhaseSpaceEventInfoSelector::operator>:" << std::endl;

//--- evaluate cut
  if ( cut_ ) {
    edm::Handle<GenPhaseSpaceEventInfo> genPhaseSpaceEventInfo;
    evt.getByLabel(src_, genPhaseSpaceEventInfo);
    if ( genPhaseSpaceEventInfo.isValid() ) {
      return (*cut_)(*genPhaseSpaceEventInfo);
    } else {
      edm::LogWarning ("GenPhaseSpaceEventInfoSelector") << " Failed to access GenPhaseSpaceEventInfo from src = " << src_.label()
							 << " --> rejecting Event !!";
      return false;
    }
  } else {
    return true;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(EventSelectorPluginFactory, GenPhaseSpaceEventInfoSelector, "GenPhaseSpaceEventInfoSelector");      
