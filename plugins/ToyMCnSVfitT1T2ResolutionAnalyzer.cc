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

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesis.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisFwd.h"

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
  edm::Handle<NSVfitEventHypothesisCollection> nSVfitEventHypotheses;
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

  for ( NSVfitEventHypothesisCollection::const_iterator nSVfitEventHypothesis = nSVfitEventHypotheses->begin();
	nSVfitEventHypothesis != nSVfitEventHypotheses->end(); ++nSVfitEventHypothesis ) {

    const edm::OwnVector<NSVfitResonanceHypothesis>& resonances = nSVfitEventHypothesis->resonances();     
    assert(resonances.size() == 1);
    if ( !resonances[0].isValidSolution() ) continue;
    const edm::OwnVector<NSVfitSingleParticleHypothesisBase>& daughters = resonances[0].daughters();
    assert(daughters.size() == 2);

    const NSVfitSingleParticleHypothesisBase* leg1 = resonances[0].daughter("leg1");
    reco::Candidate::LorentzVector p4RecLeg1 = leg1->p4();
    reco::Candidate::LorentzVector p4GenLeg1 = getGenP4<T1>(leg1);
    hGenLeg1Pt_->Fill(p4GenLeg1.pt());
    hRecLeg1Pt_->Fill(p4RecLeg1.pt());
    if ( p4GenLeg1.pt() > 0. ) hDeltaLeg1Pt_->Fill((p4RecLeg1.pt() - p4GenLeg1.pt())/p4GenLeg1.pt());

    const NSVfitSingleParticleHypothesisBase* leg2 = resonances[0].daughter("leg2");
    reco::Candidate::LorentzVector p4RecLeg2 = leg2->p4();
    reco::Candidate::LorentzVector p4GenLeg2 = getGenP4<T2>(leg2);
    hGenLeg2Pt_->Fill(p4GenLeg2.pt());
    hRecLeg2Pt_->Fill(p4RecLeg2.pt());
    if ( p4GenLeg2.pt() > 0. ) hDeltaLeg2Pt_->Fill((p4RecLeg2.pt() - p4GenLeg2.pt())/p4GenLeg2.pt());

    hDPhi_->Fill(TMath::Abs(normalizedPhi(p4RecLeg1.phi() - p4RecLeg2.phi())));

    reco::Candidate::LorentzVector p4RecMEt = nSVfitEventHypothesis->p4MEt();
    const pat::MET* toyMEt = dynamic_cast<const pat::MET*>(nSVfitEventHypothesis->met().get());
    assert(toyMEt && toyMEt->genMET());
    reco::Candidate::LorentzVector p4GenMEt = toyMEt->genMET()->p4();
    hGenMEt_->Fill(p4GenMEt.pt());
    hRecMEt_->Fill(p4RecMEt.pt());
    if ( p4GenMEt.pt() > 0. ) hDeltaMEt_->Fill((p4RecMEt.pt() - p4GenMEt.pt())/p4GenMEt.pt());
    
    hIsValidSolution_->Fill(resonances[0].isValidSolution());

    if ( resonances[0].isValidSolution() ) {
      hGenMass_->Fill(genMass);
      hRecMassMedian_->Fill(resonances[0].mass_median());
      if ( genMass > 0. ) {
	hDeltaMassMedian_->Fill((resonances[0].mass_median() - genMass)/genMass);
	hDeltaMassMedianS_->Fill((resonances[0].mass_median() - genMass)/genMass);
      }
      hRecMassMax_->Fill(resonances[0].mass_maxInterpol());
      if ( genMass > 0. ) {
	hDeltaMassMax_->Fill((resonances[0].mass_maxInterpol() - genMass)/genMass);
	hDeltaMassMaxS_->Fill((resonances[0].mass_maxInterpol() - genMass)/genMass);
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



