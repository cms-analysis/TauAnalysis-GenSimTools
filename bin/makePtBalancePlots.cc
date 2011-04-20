#include <stdio.h>
#include <iostream>
#include <sstream>
#include <memory>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include "THStack.h"

#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooLinkedList.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"

// PDFs
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooLandau.h"
#include "RooGExpModel.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooGamma.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"

using namespace RooFit;


int main(int argc, const char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: ./run_fit file_glob" << std::endl;
    return 1;
  }
  gROOT->SetBatch(true);
  // Load the data
  std::cout << "Loading data from " << argv[1] <<  std::endl;
  TChain chain("makePtBalanceNtuple/ptBalanceNtuple");
  chain.Add(argv[1]);

  TTree* dataTree = static_cast<TTree*>(&chain);

  /************************************************************************
   ********       Define variabes and data                             ****
   ************************************************************************/

  RooRealVar leg1Pt("leg1Pt", "leg1Pt", 0, 250);
  RooRealVar leg2Pt("leg2Pt", "leg2Pt", 0, 250);
  RooRealVar leg1VisPt("leg1VisPt", "leg1VisPt", 0, 200);
  RooRealVar leg2VisPt("leg2VisPt", "leg2VisPt", 0, 200);
  RooRealVar leg1VisEta("leg1VisEta", "leg1VisEta", -5, 5);
  RooRealVar leg2VisEta("leg2VisEta", "leg2VisEta", -5, 5);
  RooRealVar leg1Leg2DPhi("leg1Leg2DPhi", "leg1Leg2DPhi", 0, 3.14159);
  RooRealVar resonanceMass("resonanceMass", "resonanceMass", 0, 500);
  RooRealVar resonancePt("resonancePt", "resonancePt", 0, 500);

  TCanvas canvas("blah", "blah", 800, 600);

  RooDataSet data("data", "dataset",
      RooArgSet(
        leg1Pt, leg2Pt,
        leg1VisPt, leg2VisPt,
        leg1VisEta, leg2VisEta,
        leg1Leg2DPhi, resonanceMass, resonancePt),
      Import(*dataTree));

  RooFormulaVar scaledLeg1PtFunc("scaledLeg1Pt",
      "2.0*leg1Pt/resonanceMass", RooArgSet(leg1Pt, resonanceMass));
  RooRealVar* scaledLeg1Pt = static_cast<RooRealVar*>(
      data.addColumn(scaledLeg1PtFunc));

  std::cout << "Initial data set has " << data.numEntries() << std::endl;

  //RooAbsData* selectedData = data.reduce(
      //"leg1VisPt > 15 && leg2VisPt > 20"
      //"&& abs(leg1VisEta) < 2.1 && abs(leg2VisEta) < 2.3");
  RooAbsData* selectedData = data.reduce(
      "leg1VisPt > 0 && leg2VisPt > 0"
      "&& abs(leg1VisEta) < 2.1 && abs(leg2VisEta) < 2.3");

  std::cout << "After analysis selection data has "
    << selectedData->numEntries() << "entries" <<  std::endl;

  // Slice up the data by opening angles
  size_t nDPhiSlices = 7;
  double maxDPhi = TMath::Pi()*5.0/8;
  for(size_t i = 0; i < nDPhiSlices; i++) {
    double step = maxDPhi/nDPhiSlices;
    double maxDPhi = TMath::Pi() - (i*step);
    double minDPhi = TMath::Pi() - ((i+1)*step);
    std::stringstream dPhiCut;
    dPhiCut << "(" << minDPhi << " < leg1Leg2DPhi) && (leg1Leg2DPhi <= "
      << maxDPhi << ")";

    RooAbsData* dPhiSlice = selectedData->reduce(dPhiCut.str().c_str());
    std::cout << " Looking at a DPhi slice in (" << minDPhi << ", "
      << maxDPhi << ") with " << dPhiSlice->numEntries() << " entries. "
      << std::endl;

    std::vector<RooAbsData*> massPoints;
    std::vector<TH1*> massHistos;

    double massCutStart = 85.0;
    double massCutEnd = 305.0;
    double massCutStep = 10;
    // Slice up by resonance mass
    size_t massCounter = 0;
    for (double lowMass = massCutStart; lowMass < massCutEnd;
        lowMass += massCutStep, massCounter++ ) {
      double highMass = lowMass + massCutStep;
      std::stringstream cut;
      cut << "(" << lowMass << " < resonanceMass) && (resonanceMass <= "
        << highMass << ")";
      RooAbsData* massSlice = dPhiSlice->reduce(cut.str().c_str());
      std::cout << " --> making mass slice between (" << lowMass
        << ", " << highMass << ")" << massSlice->numEntries()
        << " entries. " << std::endl;
      std::stringstream histoname;
      histoname << "histo_mass_" << lowMass << "_" << highMass;
      TH1* histo = massSlice->createHistogram(
          histoname.str().c_str(),
          *scaledLeg1Pt,  Binning(25, 0, 5));
      histo->SetLineColor(massCounter+1);
      massPoints.push_back(massSlice);
      massHistos.push_back(histo);
    }
    THStack stack("stack", "stack");
    for (size_t iMass = 0; iMass < massPoints.size(); iMass++) {
      massHistos[iMass]->Scale(1.0/massHistos[iMass]->Integral());
      stack.Add(massHistos[iMass], "hist");
    }
    stack.Draw("nostack");

    std::stringstream file_name;
    file_name << "slice_" << i << ".pdf";
    canvas.SaveAs(file_name.str().c_str());

    for (size_t j = 0; j < massPoints.size(); j++) {
      delete massPoints[j];
      delete massHistos[j];
    }

    RooRealVar landauMean("landauMean", "landauMean", 1, 0, 10);
    RooRealVar landauWidth("landauWidth", "landauWidth", 1, 0, 10);
    RooGaussian landau("landau", "landau", *scaledLeg1Pt,
        landauMean, landauWidth);

    landau.fitTo(*dPhiSlice);
    RooPlot* scaledLeg1Frame = scaledLeg1Pt->frame(Range(0, 3.0));
    dPhiSlice->plotOn(scaledLeg1Frame);
    landau.plotOn(scaledLeg1Frame);
    scaledLeg1Frame->Draw();
    std::stringstream fit_file_name;
    fit_file_name << "slice_fit_" << i << ".pdf";
    canvas.SaveAs(fit_file_name.str().c_str());
  }
  return 0;
}
