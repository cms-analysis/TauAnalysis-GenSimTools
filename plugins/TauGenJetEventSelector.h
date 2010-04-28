//
// $Id: TauGenJetEventSelector.h,v 1.1 2009/01/21 17:32:25 veelken Exp $
//

#ifndef TauAnalysis_GenSimTools_TauGenJetEventSelector_h
#define TauAnalysis_GenSimTools_TauGenJetEventSelector_h

#include "CommonTools/UtilAlgos/interface/AnySelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectCountEventSelector.h"
#include "CommonTools/UtilAlgos/interface/MinNumberSelector.h"
#include "PhysicsTools/PatUtils/interface/MaxNumberSelector.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

typedef ObjectCountEventSelector<reco::GenJetCollection, AnySelector, MinNumberSelector> TauGenJetMinEventSelector;

typedef ObjectCountEventSelector<reco::GenJetCollection, AnySelector, MaxNumberSelector> TauGenJetMaxEventSelector;

#endif
