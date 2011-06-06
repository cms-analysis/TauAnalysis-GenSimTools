#include "TauAnalysis/GenSimTools/plugins/ToyMCnSVfitT1T2ResolutionAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h" 
#include "DataFormats/PatCandidates/interface/MET.h" 
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Math/interface/normalizedPhi.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisBase.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisBaseFwd.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesis.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesis.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesisByIntegration.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitSingleParticleHypothesis.h"

#include <TMath.h>

const double epsilon = 0.01;

template<typename T1, typename T2>
ToyMCnSVfitT1T2ResolutionAnalyzer<T1,T2>::ToyMCnSVfitT1T2ResolutionAnalyzer(const edm::ParameterSet& cfg)
{
  srcNSVfitEventHypothesis_ = cfg.getParameter<edm::InputTag>("srcNSVfitEventHypothesis");
  srcGenParticles_ = cfg.getParameter<edm::InputTag>("srcGenParticles");

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

template<typename T1, typename T2>
ToyMCnSVfitT1T2ResolutionAnalyzer<T1,T2>::~ToyMCnSVfitT1T2ResolutionAnalyzer()
{
// nothing to be done yet...
}

template<typename T1, typename T2>
void ToyMCnSVfitT1T2ResolutionAnalyzer<T1,T2>::beginJob()
{
  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_);
  
  hGenLeg1Pt_        = dqmStore.book1D("GenLeg1Pt",        "GenLeg1Pt",         50,  0.,   250.);
  hRecLeg1Pt_        = dqmStore.book1D("RecLeg1Pt",        "RecLeg1Pt",         50,  0.,   250.);
  hDeltaLeg1Pt_      = dqmStore.book1D("DeltaLeg1Pt",      "DeltaLeg1Pt",      350, -1.01,   2.49);

  hGenLeg2Pt_        = dqmStore.book1D("GenLeg2Pt",        "GenLeg2Pt",         50,  0.,   250.);
  hRecLeg2Pt_        = dqmStore.book1D("RecLeg2Pt",        "RecLeg2Pt",         50,  0.,   250.);
  hDeltaLeg2Pt_      = dqmStore.book1D("DeltaLeg2Pt",      "DeltaLeg2Pt",      350, -1.01,   2.49);
  
  hDPhi_             = dqmStore.book1D("DPhi",             "DPhi",              36, -epsilon, TMath::Pi() + epsilon);

  hGenMEt_           = dqmStore.book1D("GenMEt",           "GenMEt",            50,  0.,   250.);
  hRecMEt_           = dqmStore.book1D("RecMEt",           "RecMEt",            50,  0.,   250.);
  hDeltaMEt_         = dqmStore.book1D("DeltaMEt",         "DeltaMEt",         350, -1.01,   2.49);

  hIsValidSolution_  = dqmStore.book1D("IsValidSolution",  "IsValidSolution",    2, -0.5,    1.5);

  hGenMass_          = dqmStore.book1D("GenMass",          "GenMass",          100,  0.,   500.);
  hRecMassMedian_    = dqmStore.book1D("RecMassMedian",    "RecMassMedian",    100,  0.,   500.);
  hDeltaMassMedian_  = dqmStore.book1D("DeltaMassMedian",  "DeltaMassMedian",  350, -1.01,   2.49);
  hDeltaMassMedianS_ = dqmStore.book1D("DeltaMassMedianS", "DeltaMassMedianS",  80, -0.40,  +0.40);
  hRecMassMax_       = dqmStore.book1D("RecMassMax",       "RecMassMax",       100,  0.,   500.);
  hDeltaMassMax_     = dqmStore.book1D("DeltaMassMax",     "DeltaMassMax",     350, -1.01,   2.49);
  hDeltaMassMaxS_    = dqmStore.book1D("DeltaMassMaxS",    "DeltaMassMaxS",     80, -0.40,  +0.40);
}

//
//-------------------------------------------------------------------------------
//

template<typename T>
reco::Candidate::LorentzVector getGenP4(const NSVfitSingleParticleHypothesisBase* hypothesis)
{
  assert(0);
}

template <>
reco::Candidate::LorentzVector getGenP4<pat::Electron>(const NSVfitSingleParticleHypothesisBase* hypothesis)
{
  const pat::Electron* patElectron = dynamic_cast<const pat::Electron*>(hypothesis->particle().get());
  assert(patElectron && patElectron->genLepton());
  return patElectron->genLepton()->p4();
}

template <>
reco::Candidate::LorentzVector getGenP4<pat::Muon>(const NSVfitSingleParticleHypothesisBase* hypothesis)
{
  const pat::Muon* patMuon = dynamic_cast<const pat::Muon*>(hypothesis->particle().get());
  assert(patMuon && patMuon->genLepton());
  return patMuon->genLepton()->p4();
}

template <>
reco::Candidate::LorentzVector getGenP4<pat::Tau>(const NSVfitSingleParticleHypothesisBase* hypothesis)
{
  const pat::Tau* patTau = dynamic_cast<const pat::Tau*>(hypothesis->particle().get());
  assert(patTau && patTau->genJet());
  return patTau->genJet()->p4();
}

//
//-------------------------------------------------------------------------------
//

template<typename T1, typename T2>
void ToyMCnSVfitT1T2ResolutionAnalyzer<T1,T2>::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<NSVfitEventHypothesisBaseCollection> nSVfitEventHypotheses;
  evt.getByLabel(srcNSVfitEventHypothesis_, nSVfitEventHypotheses);

  if ( nSVfitEventHypotheses->size() != 1 ) {
    if ( nSVfitEventHypotheses->size() > 1 ) 
      edm::LogWarning("ToyMCnSVfitT1T2ResolutionAnalyzer::analyze")
	<< " Failed to find unique NSVfitEventHypothesis object !!";
    return;
  }
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByLabel(srcGenParticles_, genParticles);

  if ( genParticles->size() != 2 ) {
    edm::LogWarning("ToyMCnSVfitT1T2ResolutionAnalyzer::analyze")
      << " Failed to match NSVfitEventHypothesis object to generator level tau lepton pair !!";
    return;
  }

  double genMass = (genParticles->at(0).p4() + genParticles->at(1).p4()).mass();

  for ( NSVfitEventHypothesisBaseCollection::const_iterator nSVfitEventHypothesis = nSVfitEventHypotheses->begin();
	nSVfitEventHypothesis != nSVfitEventHypotheses->end(); ++nSVfitEventHypothesis ) {

    size_t numResonances = nSVfitEventHypothesis->numResonances();   
    assert(numResonances == 1);
    const NSVfitResonanceHypothesisBase* resonance = nSVfitEventHypothesis->resonance(0);

    hIsValidSolution_->Fill(resonance->isValidSolution());

    if ( !resonance->isValidSolution() ) continue;

    size_t numDaughters = resonance->numDaughters();   
    assert(numDaughters == 2);
    const NSVfitSingleParticleHypothesisBase* daughter1 = resonance->daughter(0);
    const NSVfitSingleParticleHypothesisBase* daughter2 = resonance->daughter(1);

    reco::Candidate::LorentzVector p4VisRecLeg1 = daughter1->particle()->p4();
    reco::Candidate::LorentzVector p4VisGenLeg1 = getGenP4<T1>(daughter1);
    hGenLeg1Pt_->Fill(p4VisGenLeg1.pt());
    hRecLeg1Pt_->Fill(p4VisRecLeg1.pt());
    if ( p4VisGenLeg1.pt() > 0. ) hDeltaLeg1Pt_->Fill((p4VisRecLeg1.pt() - p4VisGenLeg1.pt())/p4VisGenLeg1.pt());
    
    reco::Candidate::LorentzVector p4VisRecLeg2 = daughter2->particle()->p4();
    reco::Candidate::LorentzVector p4VisGenLeg2 = getGenP4<T2>(daughter2);
    hGenLeg2Pt_->Fill(p4VisGenLeg2.pt());
    hRecLeg2Pt_->Fill(p4VisRecLeg2.pt());
    if ( p4VisGenLeg2.pt() > 0. ) hDeltaLeg2Pt_->Fill((p4VisRecLeg2.pt() - p4VisGenLeg2.pt())/p4VisGenLeg2.pt());

    hDPhi_->Fill(TMath::Abs(normalizedPhi(p4VisRecLeg1.phi() - p4VisRecLeg2.phi())));
    
    reco::Candidate::LorentzVector p4RecMEt = nSVfitEventHypothesis->met()->p4();
    const pat::MET* toyMEt = dynamic_cast<const pat::MET*>(nSVfitEventHypothesis->met().get());
    assert(toyMEt && toyMEt->genMET());
    reco::Candidate::LorentzVector p4GenMEt = toyMEt->genMET()->p4();
    hGenMEt_->Fill(p4GenMEt.pt());
    hRecMEt_->Fill(p4RecMEt.pt());
    if ( p4GenMEt.pt() > 0. ) hDeltaMEt_->Fill((p4RecMEt.pt() - p4GenMEt.pt())/p4GenMEt.pt());

    if ( dynamic_cast<const NSVfitResonanceHypothesisByIntegration*>(resonance) ) {
      const NSVfitResonanceHypothesisByIntegration* resonance_byIntegration = 
	dynamic_cast<const NSVfitResonanceHypothesisByIntegration*>(resonance);
      hGenMass_->Fill(genMass);
      hRecMassMedian_->Fill(resonance_byIntegration->mass_median());
      if ( genMass > 0. ) {
	hDeltaMassMedian_->Fill((resonance_byIntegration->mass_median() - genMass)/genMass);
	hDeltaMassMedianS_->Fill((resonance_byIntegration->mass_median() - genMass)/genMass);
      }
      hRecMassMax_->Fill(resonance_byIntegration->mass_maxInterpol());
      if ( genMass > 0. ) {
	hDeltaMassMax_->Fill((resonance_byIntegration->mass_maxInterpol() - genMass)/genMass);
	hDeltaMassMaxS_->Fill((resonance_byIntegration->mass_maxInterpol() - genMass)/genMass);
      }
    }
  }
}

typedef ToyMCnSVfitT1T2ResolutionAnalyzer<pat::Electron, pat::Muon> ToyMCnSVfitElecMuPairResolutionAnalyzer;
typedef ToyMCnSVfitT1T2ResolutionAnalyzer<pat::Muon, pat::Tau> ToyMCnSVfitMuTauPairResolutionAnalyzer;
typedef ToyMCnSVfitT1T2ResolutionAnalyzer<pat::Tau, pat::Tau>ToyMCnSVfitDiTauPairResolutionAnalyzer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ToyMCnSVfitElecMuPairResolutionAnalyzer);
DEFINE_FWK_MODULE(ToyMCnSVfitMuTauPairResolutionAnalyzer);
DEFINE_FWK_MODULE(ToyMCnSVfitDiTauPairResolutionAnalyzer);



