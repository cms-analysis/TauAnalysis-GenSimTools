#include <boost/foreach.hpp>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "TTree.h"
#include "TMath.h"
#include <Math/VectorUtil.h>

#include <map>


class NSVfitPtBalanceNtupleProducer : public edm::EDAnalyzer {
  public:
    NSVfitPtBalanceNtupleProducer(const edm::ParameterSet& pset);
    virtual ~NSVfitPtBalanceNtupleProducer();
    void analyze(const edm::Event& evt, const edm::EventSetup& es);
  private:
    void addBranch(const std::string& name);
    void setValue(const std::string& name, double value);
    void fillKinematics(const std::string& prefix,
        const reco::Candidate::LorentzVector& p4);

    edm::InputTag leg1Src_;
    edm::InputTag leg2Src_;
    edm::InputTag leg1VisSrc_;
    edm::InputTag leg2VisSrc_;
    double resonanceMass_;
    TTree * ntuple_;
    typedef std::map<std::string, Double_t*> BranchMap;
    BranchMap branches_;
};

void NSVfitPtBalanceNtupleProducer::addBranch(const std::string& name) {
  std::string formatName = name + "/D";
  Double_t* destination = new Double_t;
  ntuple_->Branch(name.c_str(), destination, formatName.c_str());
  assert(branches_.count(name) == 0);
  branches_[name] = destination;
}

void NSVfitPtBalanceNtupleProducer::setValue(
    const std::string& name, double value) {
  BranchMap::iterator branchIter = branches_.find(name);
  if (branchIter != branches_.end()) {
    *(branchIter->second) = value;
  } else {
    throw cms::Exception("BadBranch") <<
      "The branch " << name << " d.n.e";
  }
}

void NSVfitPtBalanceNtupleProducer::fillKinematics(const std::string& prefix,
    const reco::Candidate::LorentzVector& p4) {
  setValue(prefix + "Pt", p4.pt());
  setValue(prefix + "Eta", p4.eta());
  setValue(prefix + "Phi", p4.phi());
  setValue(prefix + "Mass", p4.mass());
}

// Cleanup
NSVfitPtBalanceNtupleProducer::~NSVfitPtBalanceNtupleProducer() {
  BOOST_FOREACH(BranchMap::value_type ptr, branches_) {
    delete ptr.second;
  }
}

NSVfitPtBalanceNtupleProducer::NSVfitPtBalanceNtupleProducer(
    const edm::ParameterSet& pset) {
  leg1Src_ = pset.getParameter<edm::InputTag>("leg1Src");
  leg2Src_ = pset.getParameter<edm::InputTag>("leg2Src");
  leg1VisSrc_ = pset.getParameter<edm::InputTag>("leg1VisSrc");
  leg2VisSrc_ = pset.getParameter<edm::InputTag>("leg2VisSrc");
  resonanceMass_ = pset.getParameter<double>("resonanceMass");

  // Build our ntuple
  edm::Service<TFileService> fs;
  ntuple_ = fs->make<TTree>("ptBalanceNtuple", "ptBalanceNtuple");

  // Build the branches
  addBranch("leg1Pt");
  addBranch("leg1Eta");
  addBranch("leg1Phi");
  addBranch("leg1Mass");
  addBranch("leg1Charge");

  addBranch("leg1VisPt");
  addBranch("leg1VisEta");
  addBranch("leg1VisPhi");
  addBranch("leg1VisMass");

  addBranch("leg1InvisPt");
  addBranch("leg1InvisEta");
  addBranch("leg1InvisPhi");
  addBranch("leg1InvisMass");
  addBranch("leg1GJAngleLab");
  addBranch("leg1GJAngleRest");

  addBranch("leg2Pt");
  addBranch("leg2Eta");
  addBranch("leg2Phi");
  addBranch("leg2Mass");

  addBranch("leg2VisPt");
  addBranch("leg2VisEta");
  addBranch("leg2VisPhi");
  addBranch("leg2VisMass");
  addBranch("leg2Charge");

  addBranch("leg2InvisPt");
  addBranch("leg2InvisEta");
  addBranch("leg2InvisPhi");
  addBranch("leg2InvisMass");

  addBranch("leg2GJAngleLab");
  addBranch("leg2GJAngleRest");

  addBranch("leg1Leg2DPhi");

  addBranch("resonancePt");
  addBranch("resonanceEta");
  addBranch("resonancePhi");
  // These two quantities should be the same, minus the width of the resonance
  addBranch("resonanceMass");
  addBranch("resonanceInput");

  addBranch("resonanceVisPt");
  addBranch("resonanceVisEta");
  addBranch("resonanceVisPhi");
  addBranch("resonanceVisMass");

  addBranch("genMetPt");
  addBranch("genMetEta"); // not really a met quantity, but whatever
  addBranch("genMetPhi");
  addBranch("genMetMass"); // just so automated filler works

  addBranch("hadronicRecoil");

  addBranch("mtLeg1Met");
  addBranch("mtLeg2Met");

}

double mt(const reco::Candidate::LorentzVector& p4Cand,
    const reco::Candidate::LorentzVector& metP4) {
  return TMath::Sqrt(2*p4Cand.pt()*metP4.pt()*(1-TMath::Cos(
          reco::deltaPhi(p4Cand.phi(), metP4.phi()))));
}

double gjAngleLab(
    const reco::Candidate::LorentzVector& visP4,
    const reco::Candidate::LorentzVector& tauP4) {
  return std::abs(ROOT::Math::VectorUtil::Angle(tauP4, visP4));
}

double gjAngleRest(
    const reco::Candidate::LorentzVector& visP4,
    const reco::Candidate::LorentzVector& tauP4) {
  reco::Candidate::Vector boost = tauP4.BoostToCM();
  reco::Candidate::LorentzVector visP4Rest =
    ROOT::Math::VectorUtil::boost(visP4, boost);
  return std::abs(ROOT::Math::VectorUtil::Angle(tauP4, visP4Rest));
}

void NSVfitPtBalanceNtupleProducer::analyze(
    const edm::Event& evt, const edm::EventSetup& es) {

  edm::Handle<edm::View<reco::Candidate> > leg1;
  evt.getByLabel(leg1Src_, leg1);
  edm::Handle<edm::View<reco::Candidate> > leg2;
  evt.getByLabel(leg2Src_, leg2);

  edm::Handle<edm::View<reco::Candidate> > leg1Vis;
  evt.getByLabel(leg1VisSrc_, leg1Vis);

  edm::Handle<edm::View<reco::Candidate> > leg2Vis;
  evt.getByLabel(leg2VisSrc_, leg2Vis);

  reco::Candidate::LorentzVector leg1P4 = (*leg1)[0].p4();
  reco::Candidate::LorentzVector leg2P4 = (*leg2)[0].p4();
  reco::Candidate::LorentzVector leg1VisP4 = (*leg1Vis)[0].p4();
  reco::Candidate::LorentzVector leg2VisP4 = (*leg2Vis)[0].p4();
  reco::Candidate::LorentzVector leg1InvisP4 = leg1P4-leg1VisP4;
  reco::Candidate::LorentzVector leg2InvisP4 = leg2P4-leg2VisP4;

  setValue("leg1Charge", (*leg1)[0].charge());
  setValue("leg2Charge", (*leg2)[0].charge());

  // Fill the kinematic branches
  fillKinematics("leg1", leg1P4);
  fillKinematics("leg2", leg2P4);
  fillKinematics("leg1Vis", leg1VisP4);
  fillKinematics("leg2Vis", leg2VisP4);
  fillKinematics("leg1Invis", leg1InvisP4);
  fillKinematics("leg2Invis", leg2InvisP4);

  // Compute opening angles
  setValue("leg1GJAngleRest", gjAngleRest(leg1VisP4, leg1P4));
  setValue("leg1GJAngleLab", gjAngleLab(leg1VisP4, leg1P4));
  setValue("leg2GJAngleRest", gjAngleRest(leg2VisP4, leg2P4));
  setValue("leg2GJAngleLab", gjAngleLab(leg2VisP4, leg2P4));

  reco::Candidate::LorentzVector resonance = leg1P4 + leg2P4;
  fillKinematics("resonance", resonance);

  reco::Candidate::LorentzVector resonanceVis = leg1VisP4 + leg2VisP4;
  fillKinematics("resonanceVis", resonanceVis);

  setValue("leg1Leg2DPhi",
      std::abs(reco::deltaPhi(leg1VisP4.phi(), leg2VisP4.phi())));

  reco::Candidate::LorentzVector leg1Invis = leg1P4 - leg1VisP4;
  reco::Candidate::LorentzVector leg2Invis = leg2P4 - leg2VisP4;

  reco::Candidate::LorentzVector genMet = leg1Invis + leg2Invis;

  fillKinematics("genMet", genMet);

  setValue("hadronicRecoil", (genMet + resonanceVis).pt());

  setValue("mtLeg1Met", mt(leg1VisP4, genMet));
  setValue("mtLeg2Met", mt(leg2VisP4, genMet));

  ntuple_->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(NSVfitPtBalanceNtupleProducer);
