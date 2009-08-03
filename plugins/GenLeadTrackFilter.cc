// -*- C++ -*-
//
// Package:    GenLeadTrackFilter
// Class:      GenLeadTrackFilter
// 
/**\class GenLeadTrackFilter GenLeadTrackFilter.cc WToTauNuOfflineAnalysis/GenLeadTrackFilter/src/GenLeadTrackFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Souvik DAS
//         Created:  Mon Jul  6 16:21:48 CEST 2009
// $Id: GenLeadTrackFilter.cc,v 1.1 2009/07/14 12:28:26 sdas Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Particle.h"

using namespace edm;
using namespace reco;
using namespace math;

//
// class declaration
//

class GenLeadTrackFilter : public edm::EDFilter {
   public:
      explicit GenLeadTrackFilter(const edm::ParameterSet&);
      ~GenLeadTrackFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      InputTag hepMCProduct_label_;
      double   genLeadTrackPt_,
               genEta_;
      
      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenLeadTrackFilter::GenLeadTrackFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   hepMCProduct_label_   = iConfig.getParameter<InputTag>("HepMCProduct");
   genLeadTrackPt_ = iConfig.getParameter<double>("GenLeadTrackPt");
   genEta_         = iConfig.getParameter<double>("GenEta");

}


GenLeadTrackFilter::~GenLeadTrackFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool GenLeadTrackFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool allow=false;
  
  ESHandle<HepPDT::ParticleDataTable> pdt;
  iSetup.getData(pdt);

  Handle<HepMCProduct> genEvent;
  iEvent.getByLabel(hepMCProduct_label_, genEvent);
  if (genEvent.isValid())
  {
    float genLeadTrackPt=-100;
    for (HepMC::GenEvent::particle_const_iterator iter=(*(genEvent->GetEvent())).particles_begin();
         iter!=(*(genEvent->GetEvent())).particles_end(); 
         ++iter)
    {
      HepMC::GenParticle* theParticle=*iter;
      double pt=pow(pow(theParticle->momentum().px(),2)+pow(theParticle->momentum().py(),2), 0.5);
      double charge=pdt->particle(theParticle->pdg_id())->charge();
      if (theParticle->status()==1 &&
          charge!=0 &&
          fabs(theParticle->momentum().eta())<genEta_ &&
          pt>genLeadTrackPt)
      {
        genLeadTrackPt=pt;
      }
    }
    if (genLeadTrackPt>genLeadTrackPt_) allow=true; else allow=false;
  }
  else {std::cout<<"genEvent in not valid!"<<std::endl;}
  return allow;
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenLeadTrackFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenLeadTrackFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenLeadTrackFilter);
