#include "TauAnalysis/GenSimTools/plugins/ToyPATTauProducer.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

ToyPATTauProducer::ToyPATTauProducer(const edm::ParameterSet& cfg)
  : ToyParticleProducerBase(cfg)
{ 
  produces<reco::PFCandidateCollection>("");
  produces<pat::TauCollection>("");
}

ToyPATTauProducer::~ToyPATTauProducer()
{
// nothing to be done yet...
}

reco::PFCandidateRefVector getPFCandidateRefVector(reco::PFCandidateRefProd& pfCandidateRefs, const std::vector<unsigned>& indices)
{
  reco::PFCandidateRefVector retVal;
  for ( std::vector<unsigned>::const_iterator idx = indices.begin();
	idx != indices.end(); ++idx ){
    retVal.push_back(reco::PFCandidateRef(pfCandidateRefs, *idx));
  }
  return retVal;
}

void ToyPATTauProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  std::auto_ptr<reco::PFCandidateCollection> pfCandidates(new reco::PFCandidateCollection());
  reco::PFCandidateRefProd pfCandidateRefs = evt.getRefBeforePut<reco::PFCandidateCollection>();

  std::auto_ptr<pat::TauCollection> patTaus(new pat::TauCollection());

  edm::Handle<reco::GenJetCollection> genTauJets;
  evt.getByLabel(src_, genTauJets);

  unsigned numGenTauJets = genTauJets->size();
  unsigned idxPFCandidate = 0;
  for ( unsigned idxGenTauJet = 0; idxGenTauJet < numGenTauJets; ++idxGenTauJet ) {
    reco::GenJetRef genTauJet(genTauJets, idxGenTauJet);
    std::vector<const reco::GenParticle*> genDecayProducts = genTauJet->getGenConstituents();
    std::vector<unsigned> idxPFGammas;
    std::vector<unsigned> idxPFChargedHadrons;
    int idxLeadPFChargedHadron = -1;
    double leadPFChargedHadronPt = -1.;
    std::vector<unsigned> idxPFCandidates;
    reco::Candidate::LorentzVector p4smearedTauJet;
    double chargeTauJet = 0.;
    double rndSmearFactor = rnd_.Gaus(1., resolution_);
    for ( std::vector <const reco::GenParticle*>::const_iterator genDecayProduct = genDecayProducts.begin();
	  genDecayProduct != genDecayProducts.end(); ++genDecayProduct ) {
      int pdgId = TMath::Abs((*genDecayProduct)->pdgId());
      double charge = (*genDecayProduct)->charge();
      if ( pdgId >= 11 && pdgId <= 16 ) continue;
      reco::Candidate::LorentzVector p4smeared = (*genDecayProduct)->p4();
      p4smeared *= rndSmearFactor;
      reco::PFCandidate pfCandidate;
      pfCandidate.setP4(p4smeared);  
      if ( pdgId == 22 ) {
	pfCandidate.setParticleType(reco::PFCandidate::gamma); 
	idxPFGammas.push_back(idxPFCandidate);
      } else if ( charge != 0. ) {
	pfCandidate.setParticleType(reco::PFCandidate::h); 
	idxPFChargedHadrons.push_back(idxPFCandidate);
	if ( pfCandidate.pt() > leadPFChargedHadronPt ) {
	  idxLeadPFChargedHadron = idxPFCandidate;
	  leadPFChargedHadronPt = pfCandidate.pt();
	}
      } 
      pfCandidates->push_back(pfCandidate);
      idxPFCandidates.push_back(idxPFCandidate);
      p4smearedTauJet += p4smeared;
      chargeTauJet += charge;
      ++idxPFCandidate;
    } 
    assert(idxLeadPFChargedHadron != -1);
    assert(chargeTauJet == -1. || chargeTauJet == +1.);
    reco::PFTau tau;
    tau.setP4(p4smearedTauJet);
    tau.setCharge(chargeTauJet);
    tau.setleadPFChargedHadrCand(reco::PFCandidateRef(pfCandidateRefs, idxLeadPFChargedHadron));
    tau.setsignalPFCands(getPFCandidateRefVector(pfCandidateRefs, idxPFCandidates));
    tau.setsignalPFChargedHadrCands(getPFCandidateRefVector(pfCandidateRefs, idxPFChargedHadrons));
    tau.setsignalPFGammaCands(getPFCandidateRefVector(pfCandidateRefs, idxPFGammas));
    pat::Tau patTau(tau);
    int tauDecayMode = reco::PFTauDecayMode::tauDecayOther;
    if ( idxPFChargedHadrons.size() == 1 ) {
      if      ( idxPFGammas.size() < 2 ) tauDecayMode = reco::PFTauDecayMode::tauDecay1ChargedPion0PiZero;
      else if ( idxPFGammas.size() < 4 ) tauDecayMode = reco::PFTauDecayMode::tauDecay1ChargedPion1PiZero;
      else if ( idxPFGammas.size() < 6 ) tauDecayMode = reco::PFTauDecayMode::tauDecay1ChargedPion2PiZero;
    } else if ( idxPFChargedHadrons.size() == 3 ) {
      if      ( idxPFGammas.size() < 2 ) tauDecayMode = reco::PFTauDecayMode::tauDecay3ChargedPion0PiZero;
      else if ( idxPFGammas.size() < 4 ) tauDecayMode = reco::PFTauDecayMode::tauDecay3ChargedPion1PiZero;
    }
    patTau.setDecayMode(tauDecayMode);    
    patTau.setGenJet(genTauJet);
    //std::cout << " genTauJet: Pt = " << genTauJet->pt() << "," 
    //	        << " eta = " << genTauJet->eta() << ", phi = " << genTauJet->phi() << std::endl;
    //std::cout << "--> patTau: Pt = " << patTau.pt() << "," 
    //          << " eta = " << patTau.eta() << ", phi = " << patTau.phi() << std::endl;
    patTaus->push_back(patTau);
  }

  evt.put(pfCandidates);
  evt.put(patTaus);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ToyPATTauProducer);
