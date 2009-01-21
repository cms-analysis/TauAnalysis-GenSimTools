#include "TauAnalysis/GenSimTools/plugins/TauGenJetEventSelector.h"

#include "PhysicsTools/UtilAlgos/interface/EventSelectorBase.h"

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(EventSelectorPluginFactory, TauGenJetMinEventSelector, "TauGenJetMinEventSelector");

DEFINE_EDM_PLUGIN(EventSelectorPluginFactory, TauGenJetMaxEventSelector, "TauGenJetMaxEventSelector");
