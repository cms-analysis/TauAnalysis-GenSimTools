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
#include "RooDataHist.h"

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
#include "RooMsgService.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitPtBalancePdfs.h"

using namespace RooFit;

void plotSlices(const RooRealVar& toPlot,
    const RooAbsPdf& fittedModel, RooAbsData& data,
    const RooArgSet& conditionalVars, TPad* pad, const std::string& prefix,
    size_t massSlices, double massLow, double massHigh,
    size_t phiSlices, double phiLow, double phiHigh) {

  pad->cd();

  std::cout << "Plotting " << massSlices << " mass slices and "
    << phiSlices << " phi slices" << std::endl;

  double phiSliceSize = (phiHigh - phiLow)/phiSlices;
  double massSliceSize = (massHigh - massLow)/massSlices;
  for (size_t iMass = 0; iMass < massSlices; ++iMass) {

    double minMass = massLow + iMass*massSliceSize;
    double maxMass = massLow + (iMass+1)*massSliceSize;

    std::stringstream massSliceCut;
    massSliceCut << minMass << " < resonanceMass && resonanceMass < " << maxMass;

    std::auto_ptr<RooAbsData> massSlice(
        data.reduce(massSliceCut.str().c_str()));
    std::stringstream massRangeStr;
    massRangeStr << "(" << minMass << " - " << maxMass << ")";
    std::cout << " In mass slice " << iMass << massRangeStr.str() <<
      "there are " << massSlice->numEntries() << " entries." << std::endl;

    for (size_t iPhi = 0; iPhi < phiSlices; ++iPhi) {
      double minPhi = phiLow + iPhi*phiSliceSize;
      double maxPhi = phiLow + (iPhi+1)*phiSliceSize;
      std::stringstream phiSliceCut;
      phiSliceCut << minPhi << " < dPhi && dPhi < " << maxPhi;
      std::auto_ptr<RooAbsData> phiSlice(
          massSlice->reduce(phiSliceCut.str().c_str()));
      std::stringstream phiRangeStr;
      phiRangeStr << "(" << minPhi << " - " << maxPhi << ")";
      std::cout << " In phi slice " << iPhi << phiRangeStr.str() <<
        "there are " << phiSlice->numEntries() << " entries." << std::endl;

      std::auto_ptr<RooPlot> frame(toPlot.frame(Range(0, 4.0)));
      phiSlice->plotOn(frame.get());
      fittedModel.plotOn(frame.get(), ProjWData(conditionalVars, *phiSlice));
      std::stringstream plotTitle;
      plotTitle << " Leg1Pt for mass " << massRangeStr.str()
        << " and #Delta#phi " << phiRangeStr.str();
      frame->SetName(plotTitle.str().c_str());
      frame->Draw();
      std::stringstream sliceFileName;
      sliceFileName << prefix << "_phi_"
        << iPhi << "_mass_" << iMass << ".png";
      pad->SaveAs(sliceFileName.str().c_str());
    }
  }
}

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
  int filesAdded = chain.Add(argv[1]);

  if (!filesAdded) {
    std::cerr << "No valid input files specified!" << std::endl;
    return 1;
  }

  TTree* dataTree = static_cast<TTree*>(&chain);

  RooMsgService::instance().Print();
  std::cout << endl;
  RooMsgService::instance().addStream(DEBUG,Topic(Minimization));
  RooMsgService::instance().Print();
  std::cout << endl;

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
  scaledLeg1Pt->setRange(0, 5);
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
  std::string chebyshevFirstOrder = "abs(@1 + @2*(@0))";

  RooRealVar gaussMeanA("gaussMeanA", "gaussMeanA", 1, 0, 5);
  RooRealVar gaussMeanB("gaussMeanB", "gaussMeanB", 0, -5, 5);
  RooRealVar gaussMeanC("gaussMeanC", "gaussMeanC", 0, -5, 5);
  RooFormulaVar gaussMeanDPhi(
      "gaussMeanDPhi", "Gaussian #mu", chebyshevSecondOrder.c_str(),
      RooArgList(*dPhi, gaussMeanA, gaussMeanB, gaussMeanC));

  RooRealVar gaussSigmaA("gaussSigmaA", "gaussSigmaA", 0.1, 0, 5);
  RooRealVar gaussSigmaB("gaussSigmaB", "gaussSigmaB", 0.1, -5, 5);
  RooRealVar gaussSigmaC("gaussSigmaC", "gaussSigmaC", 0, -5, 5);
  RooFormulaVar gaussSigmaDPhi(
      "gaussSigmaDPhi", "Gaussian #sigma", chebyshevSecondOrder.c_str(),
      RooArgList(*dPhi, gaussSigmaA, gaussSigmaB, gaussSigmaC));

  RooRealVar gammaScaleA("gammaScaleA", "gammaScaleA", 10, 0, 25);
  RooRealVar gammaScaleB("gammaScaleB", "gammaScaleB", -2.5, -500, 500);
  RooRealVar gammaScaleC("gammaScaleC", "gammaScaleC", 0, -5, 5);
  RooFormulaVar gammaScaleDPhi(
      "gammaScaleDPhi", "#Gamma scale", chebyshevSecondOrder.c_str(),
      RooArgList(*dPhi, gammaScaleA, gammaScaleB, gammaScaleC));

  RooRealVar gammaShapeA("gammaShapeA", "gammaShapeA", 0.05, 0, 25);
  RooRealVar gammaShapeB("gammaShapeB", "gammaShapeB", 0.113, -25, 25);
  RooRealVar gammaShapeC("gammaShapeC", "gammaShapeC", 0, -5, 5);
  RooFormulaVar gammaShapeDPhi(
      "gammaShapeDPhi", "#Gamma scale", chebyshevSecondOrder.c_str(),
      RooArgList(*dPhi, gammaShapeA, gammaShapeB, gammaShapeC));

  RooRealVar gammaFracA("gammaFracA", "gammaFracA", 0.0096, -15, 15);
  RooRealVar gammaFracB("gammaFracB", "gammaFracB", 0.015, -15, 15);
  RooFormulaVar gammaFracDPhi(
      "gammaFracDPhi", "#Gamma frac term", chebyshevFirstOrder.c_str(),
      RooArgList(*dPhi, gammaFracA, gammaFracB));

  // Make mass based correction factors
  RooRealVar gaussMeanCorrA("gaussMeanCorrA", "gaussMeanCorrA", 1, -5, 5);
  RooRealVar gaussMeanCorrB("gaussMeanCorrB", "gaussMeanCorrB", 0, -5, 5);
  RooRealVar gaussMeanCorrC("gaussMeanCorrC", "gaussMeanCorrC", 0, -5, 5);
  RooFormulaVar gaussMeanCorr(
      "gaussMeanCorr", "Gaussian #mu", chebyshevSecondOrder.c_str(),
      RooArgList(resonanceMass, gaussMeanCorrA, gaussMeanCorrB, gaussMeanCorrC));

  RooRealVar gaussSigmaCorrA("gaussSigmaCorrA", "gaussSigmaCorrA", 1.0, -5, 5);
  RooRealVar gaussSigmaCorrB("gaussSigmaCorrB", "gaussSigmaCorrB", 0.0, -5, 5);
  RooRealVar gaussSigmaCorrC("gaussSigmaCorrC", "gaussSigmaCorrC", 0.0, -5, 5);
  RooFormulaVar gaussSigmaCorr(
      "gaussSigmaCorr", "Gaussian #sigma", chebyshevSecondOrder.c_str(),
      RooArgList(resonanceMass, gaussSigmaCorrA, gaussSigmaCorrB, gaussSigmaCorrC));

  RooRealVar gammaScaleCorrA("gammaScaleCorrA", "gammaScaleCorrA", 1.0, -25, 25);
  RooRealVar gammaScaleCorrB("gammaScaleCorrB", "gammaScaleCorrB", 0.0, -25, 25);
  RooRealVar gammaScaleCorrC("gammaScaleCorrC", "gammaScaleCorrC", 0.0, -25, 25);
  RooFormulaVar gammaScaleCorr(
      "gammaScaleCorr", "#Gamma scale", chebyshevSecondOrder.c_str(),
      RooArgList(resonanceMass, gammaScaleCorrA, gammaScaleCorrB, gammaScaleCorrC));


  RooRealVar gammaShapeCorrA("gammaShapeCorrA", "gammaShapeCorrA", 1.0, -25, 25);
  RooRealVar gammaShapeCorrB("gammaShapeCorrB", "gammaShapeCorrB", 0.0, -25, 25);
  RooRealVar gammaShapeCorrC("gammaShapeCorrC", "gammaShapeCorrC", 0.0, -25, 25);
  RooFormulaVar gammaShapeCorr(
      "gammaShapeCorr", "#Gamma scale", chebyshevSecondOrder.c_str(),
      RooArgList(resonanceMass, gammaShapeCorrA, gammaShapeCorrB, gammaShapeCorrC));

  RooRealVar gammaFracCorrA("gammaFracCorrA", "gammaFracCorrA", 1.0, -5, 5);
  RooRealVar gammaFracCorrB("gammaFracCorrB", "gammaFracCorrB", 0.0, -5, 5);
  RooFormulaVar gammaFracCorr(
      "gammaFracCorr", "#Gamma frac term", chebyshevFirstOrder.c_str(),
      RooArgList(resonanceMass, gammaFracCorrA, gammaFracCorrB));

  std::string multiplyTwo = "@0*@1";
  RooFormulaVar gaussMean(
      "gaussMean", "Gaussian #mu", multiplyTwo.c_str(),
      RooArgList(gaussMeanDPhi, gaussMeanCorr));
  RooFormulaVar gaussSigma(
      "gaussSigma", "Gaussian #sigma", multiplyTwo.c_str(),
      RooArgList(gaussSigmaDPhi, gaussSigmaCorr));
  RooFormulaVar gammaScale(
      "gammaScale", "#Gamma scale", multiplyTwo.c_str(),
      RooArgList(gammaScaleDPhi, gammaScaleCorr));
  RooFormulaVar gammaShape(
      "gammaShape", "#Gamma shape", multiplyTwo.c_str(),
      RooArgList(gammaShapeDPhi, gammaShapeCorr));

  RooFormulaVar gammaFrac(
      "gammaFrac", "#Gamma scale", "(1.0/3.14159)*atan(@0*@1)+0.5",
      RooArgList(gammaFracCorr, gammaFracDPhi));

  RooFormulaVar gammaFracForDPhi(
      "gammaFracForDPhi", "#Gamma scale", "(1.0/3.14159)*atan(@0)+0.5",
      RooArgList(gammaFracDPhi));


  // Fix mass corrections
  gaussMeanCorrA.setConstant(true);
  gaussMeanCorrB.setConstant(true);
  gaussMeanCorrC.setConstant(true);

  gaussSigmaCorrA.setConstant(true);
  gaussSigmaCorrB.setConstant(true);
  gaussSigmaCorrC.setConstant(true);

  gammaScaleCorrA.setConstant(true);
  gammaScaleCorrB.setConstant(true);
  gammaScaleCorrC.setConstant(true);

  gammaShapeCorrA.setConstant(true);
  gammaShapeCorrB.setConstant(true);
  gammaShapeCorrC.setConstant(true);

  gammaFracCorrA.setConstant(true);
  gammaFracCorrB.setConstant(true);

  RooRealVar zero("zero", "zero", 0, -1, 1);
  zero.setConstant(true);

  RooGamma gamma("gamma", "gamma",
      *scaledLeg1Pt, gammaScale, gammaShape, zero);
  RooGaussian gauss("gauss", "gauss",
      *scaledLeg1Pt, gaussMean, gaussSigma);

  // Build model
  RooAddPdf model("sum", "sum", gamma, gauss, gammaFrac);

  // Build a DPhi only model
  RooGamma gammaDPhi("gammaDPhi", "gamma",
      *scaledLeg1Pt, gammaScaleDPhi, gammaShapeDPhi, zero);

  RooGaussian gaussDPhi("gaussDPhi", "gauss",
      *scaledLeg1Pt, gaussMeanDPhi, gaussSigmaDPhi);

  // Build model
  RooAddPdf modelDPhi("sumDPhi", "sumDPhi", gammaDPhi, gaussDPhi, gammaFracForDPhi);

  model.Print("v");

  modelDPhi.Print("v");

  std::cout << "Fitting reduced model..." << std::endl;
  RooArgSet conditionalVariables(*dPhi, resonanceMass);
  std::cout << "Begining fit" << std::endl;
  modelDPhi.fitTo(*reducedDataSet, ConditionalObservables(*dPhi),
      NumCPU(6));

  RooDataHist binnedReducedData("binnedReducedData", "binnedReducedData",
      RooArgSet(*dPhi, *scaledLeg1Pt, resonanceMass), *reducedDataSet);

  //plotSlices(*scaledLeg1Pt, modelDPhi, binnedReducedData, conditionalVariables,
  //    &canvas, "reduced_fit_reduced_data", 1, 155, 165, 3, 0, 3.14/2);

  std::cout << "Creating binned dataset" << std::endl;

  RooDataHist binnedData("binnedData", "binnedData",
      RooArgSet(*dPhi, *scaledLeg1Pt, resonanceMass), *selectedData);

  //plotSlices(*scaledLeg1Pt, modelDPhi, binnedData, conditionalVariables,
  //    &canvas, "reduced_fit", 3, 100, 300, 3, 0, 3.14/2);

  std::cout << "Now fitting full model...." << std::endl;
  // Fix mass corrections
  //gaussMeanCorrA.setConstant(false);
  gaussMeanCorrB.setConstant(false);
  gaussMeanCorrB.setError(0.001);
  //gaussMeanCorrC.setConstant(false);

  //gaussSigmaCorrA.setConstant(false);
  gaussSigmaCorrB.setConstant(false);
  gaussSigmaCorrB.setError(0.001);
  //gaussSigmaCorrC.setConstant(false);

  //gammaScaleCorrA.setConstant(false);
  gammaScaleCorrB.setConstant(false);
  gammaScaleCorrB.setError(0.001);
  //gammaScaleCorrC.setConstant(false);

  //gammaShapeCorrA.setConstant(false);
  gammaShapeCorrB.setConstant(false);
  gammaShapeCorrB.setError(0.001);
  //gammaShapeCorrC.setConstant(false);

  //gammaFracCorrA.setConstant(false);
  gammaFracCorrB.setConstant(false);
  gammaFracCorrB.setError(0.001);

  model.fitTo(binnedData, ConditionalObservables(conditionalVariables),
      NumCPU(5));

  plotSlices(*scaledLeg1Pt, model, binnedData, conditionalVariables,
      &canvas, "full_fit", 5, 100, 300, 5, 0, 3.14/2);

  return 0;
}
