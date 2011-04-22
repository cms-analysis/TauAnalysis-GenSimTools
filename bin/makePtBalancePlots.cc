#include <stdio.h>
#include <iostream>
#include <sstream>
#include <memory>
#include <boost/foreach.hpp>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include "THStack.h"
#include "TGraphErrors.h"

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

#include "TauAnalysis/CandidateTools/interface/NSVfitPtBalancePdfs.h"

using namespace RooFit;

int main(int argc, const char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: ./run_fit file_glob" << std::endl;
    return 1;
  }
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
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
  RooFormulaVar dPhiFormula("dPhi", "#pi - #Delta #phi",
      "3.14159265 - leg1Leg2DPhi", RooArgSet(leg1Leg2DPhi));

  RooRealVar* scaledLeg1Pt = static_cast<RooRealVar*>(
      data.addColumn(scaledLeg1PtFunc));
  RooRealVar* dPhi = static_cast<RooRealVar*>(
      data.addColumn(dPhiFormula));
  dPhi->setRange(0, 1.8);

  std::cout << "Initial data set has " << data.numEntries() << std::endl;

  // Apply analysys selections and reduce the total # of variables
  RooAbsData* selectedData = data.reduce(
      RooArgSet(*scaledLeg1Pt, *dPhi, resonanceMass),
      "leg1VisPt > 15 && leg2VisPt > 20"
      "&& abs(leg1VisEta) < 2.1 && abs(leg2VisEta) < 2.3");

  std::cout << "After analysis selection data has "
    << selectedData->numEntries() << "entries" <<  std::endl;

  RooAbsData* reducedDataSet = selectedData->reduce(
      "165 > resonanceMass && resonanceMass > 155");

  std::cout << "Making reduced data set with"
    << reducedDataSet->numEntries() << "entries" <<  std::endl;

  std::string chebyshevSecondOrder = "abs(@1 + @2*(@0) + @3*(2*@0*@0 - 1))";

  RooRealVar gaussMeanA("gaussMeanA", "gaussMeanA", 1, 0, 5);
  RooRealVar gaussMeanB("gaussMeanB", "gaussMeanB", 0, -5, 5);
  RooRealVar gaussMeanC("gaussMeanC", "gaussMeanC", 0, -5, 5);
  RooFormulaVar gaussMean(
      "gaussMean", "Gaussian #mu", chebyshevSecondOrder.c_str(),
      RooArgList(*dPhi, gaussMeanA, gaussMeanB, gaussMeanC));

  RooRealVar gaussSigmaA("gaussSigmaA", "gaussSigmaA", 0.1, 0, 5);
  RooRealVar gaussSigmaB("gaussSigmaB", "gaussSigmaB", 0.1, -5, 5);
  RooRealVar gaussSigmaC("gaussSigmaC", "gaussSigmaC", 0, -5, 5);
  RooFormulaVar gaussSigma(
      "gaussSigma", "Gaussian #sigma", chebyshevSecondOrder.c_str(),
      RooArgList(*dPhi, gaussSigmaA, gaussSigmaB, gaussSigmaC));

  RooRealVar gammaScaleA("gammaScaleA", "gammaScaleA", 10, 0, 25);
  RooRealVar gammaScaleB("gammaScaleB", "gammaScaleB", -2.5, -25, 25);
  RooRealVar gammaScaleC("gammaScaleC", "gammaScaleC", 0, -5, 5);
  RooFormulaVar gammaScale(
      "gammaScale", "#Gamma scale", chebyshevSecondOrder.c_str(),
      RooArgList(*dPhi, gammaScaleA, gammaScaleB, gammaScaleC));

  RooRealVar gammaShapeA("gammaShapeA", "gammaShapeA", 0.05, 0, 25);
  RooRealVar gammaShapeB("gammaShapeB", "gammaShapeB", 0.113, -25, 25);
  RooRealVar gammaShapeC("gammaShapeC", "gammaShapeC", 0, -5, 5);
  RooFormulaVar gammaShape(
      "gammaShape", "#Gamma scale", chebyshevSecondOrder.c_str(),
      RooArgList(*dPhi, gammaShapeA, gammaShapeB, gammaShapeC));

  RooRealVar gammaFracA("gammaFracA", "gammaFracA", 0.0096, 0, 5);
  RooRealVar gammaFracB("gammaFracB", "gammaFracB", 0.015, 0, 5);
  //RooRealVar gammaFracC("gammaFracC", "gammaFracC", 0, -5, 5);
  RooFormulaVar gammaFrac(
      "gammaFrac", "#Gamma scale", "(2.0/3.14159)*atan(@2*@0 + @1)",
      RooArgList(*dPhi, gammaFracA, gammaFracB));

  RooRealVar zero("zero", "zero", 0, -1, 1);
  zero.setConstant(true);

  RooGamma gamma("gamma", "gamma",
      *scaledLeg1Pt, gammaScale, gammaShape, zero);

  RooGaussian gauss("gauss", "gauss",
      *scaledLeg1Pt, gaussMean, gaussSigma);

  RooAddPdf model("sum", "sum", gamma, gauss, gammaFrac);

  model.Print("v");

  std::cout << "Fitting reduced model..." << std::endl;
  model.fitTo(*reducedDataSet, ConditionalObservables(*dPhi));

  std::cout << "Now fitting full model...." << std::endl;
  model.fitTo(*selectedData, ConditionalObservables(*dPhi));

  RooPlot* scaledLeg1Frame = scaledLeg1Pt->frame(Range(0, 3.0));
  std::cout << "Plotting data..." << std::endl;
  selectedData->plotOn(scaledLeg1Frame);
  std::cout << "Plotting model..." << std::endl;
  model.plotOn(scaledLeg1Frame, ProjWData(*dPhi, *selectedData));
  //model.plotOn(scaledLeg1Frame, ProjWData(*dPhi, *selectedData),
      //Components(gauss), LineStyle(kDashed));
  scaledLeg1Frame->Draw();
  canvas.SaveAs("overall_fit.png");

  delete scaledLeg1Frame;

  int nMassSlices = 5;
  double lowMass = 90;
  double highMass = 300;
  for (int iMassSlice = 0; iMassSlice < nMassSlices; ++iMassSlice) {
    double minMass = lowMass + iMassSlice*(highMass - lowMass)/nMassSlices;
    double maxMass = lowMass + (iMassSlice+1)*(highMass - lowMass)/nMassSlices;
    std::stringstream massSliceCut;
    massSliceCut << minMass << " < resonanceMass && resonanceMass < " << maxMass;
    std::auto_ptr<RooAbsData> massSlice(
        selectedData->reduce(massSliceCut.str().c_str()));
    std::cout << "Made mass slice with " << massSlice->numEntries()
      << " entries" << std::endl;

    int nSlices = 4;
    double phiSliceSize = TMath::Pi()/10.0;
    for (int iSlice = 0; iSlice < nSlices; ++iSlice) {
      std::cout << "Plotting slice " << iSlice << std::endl;
      double minPhi = (iSlice)*phiSliceSize;
      double maxPhi = (iSlice+1)*phiSliceSize;
      std::auto_ptr<RooPlot> frame(scaledLeg1Pt->frame(Range(0, 3.0)));
      std::stringstream sliceCut;
      sliceCut << minPhi << " < dPhi && dPhi < " << maxPhi;
      std::auto_ptr<RooAbsData> dataSlice(
          massSlice->reduce(sliceCut.str().c_str()));
      dataSlice->plotOn(frame.get());
      model.plotOn(frame.get(), ProjWData(*dPhi, *dataSlice));
      frame->Draw();
      std::stringstream sliceFileName;
      sliceFileName << "slice_after_fit_" << iSlice << "_mass_" << iMassSlice << ".png";
      canvas.SaveAs(sliceFileName.str().c_str());
    }
  }

  // Add corrections based on resonance mass
  std::string doubleChebyshevSecondOrder =
    "abs(@1 + @2*(@0) + @3*(2*@0*@0 - 1) + @4*(@5) + @6*(2*@4*@4-1))";

  RooRealVar gaussMeanCorrA("gaussMeanCorrA", "gaussMeanCorrA", 0, -5, 5);
  RooRealVar gaussMeanCorrB("gaussMeanCorrB", "gaussMeanCorrB", 0, -5, 5);
  RooFormulaVar gaussMeanCorr(
      "gaussMeanCorr", "Gaussian #mu", doubleChebyshevSecondOrder.c_str(),
      RooArgList(*dPhi, gaussMeanA, gaussMeanB, gaussMeanC,
        resonanceMass, gaussMeanCorrA, gaussMeanCorrB));

  RooRealVar gaussSigmaCorrA("gaussSigmaCorrA", "gaussSigmaCorrA", 0.0, -5, 5);
  RooRealVar gaussSigmaCorrB("gaussSigmaCorrB", "gaussSigmaCorrB", 0.0, -5, 5);
  RooFormulaVar gaussSigmaCorr(
      "gaussSigmaCorr", "Gaussian #sigma", doubleChebyshevSecondOrder.c_str(),
      RooArgList(*dPhi, gaussSigmaA, gaussSigmaB, gaussSigmaC,
        resonanceMass, gaussSigmaCorrA, gaussSigmaCorrB));

  RooRealVar gammaScaleCorrA("gammaScaleCorrA", "gammaScaleCorrA", 0.0, -25, 25);
  RooRealVar gammaScaleCorrB("gammaScaleCorrB", "gammaScaleCorrB", 0.0, -25, 25);
  RooFormulaVar gammaScaleCorr(
      "gammaScaleCorr", "#Gamma scale", doubleChebyshevSecondOrder.c_str(),
      RooArgList(*dPhi, gammaScaleA, gammaScaleB, gammaScaleC,
        resonanceMass, gammaScaleCorrA, gammaScaleCorrB));

  RooRealVar gammaShapeCorrA("gammaShapeCorrA", "gammaShapeCorrA", 0.0, -25, 25);
  RooRealVar gammaShapeCorrB("gammaShapeCorrB", "gammaShapeCorrB", 0.0, -25, 25);
  RooFormulaVar gammaShapeCorr(
      "gammaShapeCorr", "#Gamma scale", doubleChebyshevSecondOrder.c_str(),
      RooArgList(*dPhi, gammaShapeA, gammaShapeB, gammaShapeC,
        resonanceMass, gammaShapeCorrA, gammaShapeCorrB));

  RooGamma correctedGamma("correctedGamma", "correctedGamma",
      *scaledLeg1Pt, gammaScaleCorr, gammaShapeCorr, zero);

  RooGaussian correctedGaussian("correctedGaussian", "correctedGaussian",
      *scaledLeg1Pt, gaussMeanCorr, gaussSigmaCorr);

  RooAddPdf correctedModel("sumCorr", "sumCorr", correctedGamma,
      correctedGaussian, gammaFrac);

  correctedModel.Print("v");

  std::cout << "Fitting corrected data..." << std::endl;
  correctedModel.fitTo(*selectedData,
      ConditionalObservables(RooArgSet(*dPhi, resonanceMass)));

  scaledLeg1Frame = scaledLeg1Pt->frame(Range(0, 3.0));
  std::cout << "Plotting corrected data..." << std::endl;
  selectedData->plotOn(scaledLeg1Frame);
  std::cout << "Plotting model..." << std::endl;
  correctedModel.plotOn(scaledLeg1Frame, ProjWData(
        RooArgSet(*dPhi, resonanceMass), *selectedData));
  scaledLeg1Frame->Draw();
  canvas.SaveAs("fitCorrected.png");

  return 0;
}
