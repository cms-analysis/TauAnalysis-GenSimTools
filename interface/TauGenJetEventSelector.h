//
// $Id: TauGenJetEventSelector.h,v 1.1 2008/03/06 09:23:10 veelken Exp $
//

#ifndef TauAnalysis_GenSimTools_TauGenJetEventSelector_h
#define TauAnalysis_GenSimTools_TauGenJetEventSelector_h

#include "PhysicsTools/UtilAlgos/interface/AnySelector.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectCountEventSelector.h"
#include "PhysicsTools/UtilAlgos/interface/MinNumberSelector.h"
#include "PhysicsTools/PatUtils/interface/MaxNumberSelector.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

typedef ObjectCountEventSelector<reco::GenJetCollection, AnySelector, MinNumberSelector> TauGenJetMinEventSelector;

typedef ObjectCountEventSelector<reco::GenJetCollection, AnySelector, MaxNumberSelector> TauGenJetMaxEventSelector;

#endif
