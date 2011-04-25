#include <stdio.h>
#include <iostream>
#include <cmath>
#include <sstream>
#include <memory>
#include <boost/foreach.hpp>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"

#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooLinkedList.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"

// PDFs
#include "RooGaussian.h"
#include "RooBifurGauss.h"
#include "RooGamma.h"
#include "RooMsgService.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitPtBalancePdfs.h"

using namespace RooFit;

class RooDoublePolyVar : public RooAbsReal {
  public:
    RooDoublePolyVar(const char* name, const char* title,
        RooAbsReal& x, RooAbsReal& y, RooArgList& coeffs,
        size_t orderX, size_t orderY, bool takeAbs):
      RooAbsReal(name, title),
      x_("x", "Dependent x", this, x),
      y_("y", "Dependent y", this, y),
      coeffs_("coeffs", "List of coefficients", this),
      orderX_(orderX), orderY_(orderY),takeAbs_(takeAbs) {
        std::cout << "RooDoublePolyVar::ctor"
          << " order X: " << orderX_
          << " order Y: " << orderY_
          << " Ncoeffecients: " << coeffs.getSize() << std::endl;
        assert((orderX_+1)*(orderY_+1) == (size_t)coeffs.getSize());
        for (int i = 0; i < coeffs.getSize(); i++) {
          RooAbsArg* toAdd = coeffs.at(i);
          std::cout << "Adding coefficient: " << toAdd->GetName() << std::endl;
          coeffs_.add(*toAdd);
        }
      }

    // Copy ctor
    RooDoublePolyVar(const RooDoublePolyVar& other,
        const char* name = 0) : RooAbsReal(other, name),
      x_("x", this, other.x_),
      y_("y", this, other.y_),
      coeffs_("coeffs", this, other.coeffs_),
      orderX_(other.orderX_),
      orderY_(other.orderY_),
      takeAbs_(other.takeAbs_) {}

    TObject* clone(const char* newName) const {
      return new RooDoublePolyVar(*this, newName);
    }

    Double_t evaluate() const {
      Double_t sum = 0;
      Double_t x = x_;
      Double_t y = y_;
      size_t coeffIndex = 0;
      for (size_t iX = 0; iX < (orderX_+1); ++iX) {
        for (size_t iY = 0; iY < (orderY_+1); ++iY) {
          const RooRealVar* coeff = dynamic_cast<const RooRealVar*>(
              coeffs_.at(coeffIndex));
          assert(coeff);
          sum += coeff->getVal()*std::pow(x, iX)*std::pow(y, iY);
          coeffIndex++;
        }
      }
      return takeAbs_ ? std::abs(sum) : sum;
    }

    // Get all the terms that are effected by X
    RooArgList termsProportionalToX() const {
      RooArgList output;
      size_t coeffIndex = 0;
      for (size_t iX = 0; iX < (orderX_+1); ++iX) {
        for (size_t iY = 0; iY < (orderY_+1); ++iY) {
          RooRealVar* coeff = dynamic_cast<RooRealVar*>(
              coeffs_.at(coeffIndex));
          assert(coeff);
          coeffIndex++;
          // Only add to the output list if there is no Y dependence
          if (iX != 0) {
            output.add(*coeff);
          }
        }
      }
      return output;
    }

    RooArgList termsProportionalToY() const {
      RooArgList output;
      size_t coeffIndex = 0;
      for (size_t iX = 0; iX < (orderX_+1); ++iX) {
        for (size_t iY = 0; iY < (orderY_+1); ++iY) {
          RooRealVar* coeff = dynamic_cast<RooRealVar*>(
              coeffs_.at(coeffIndex));
          assert(coeff);
          coeffIndex++;
          // Only add to the output list if there is no Y dependence
          if (iY != 0) {
            output.add(*coeff);
          }
        }
      }
      return output;
    }

    RooRealVar& getTerm(size_t i) const {
      RooRealVar* constant = dynamic_cast<RooRealVar*>(coeffs_.at(i));
      assert(constant);
      return *constant;
    }

  private:
    RooRealProxy x_;
    RooRealProxy y_;
    RooListProxy coeffs_;
    size_t orderX_;
    size_t orderY_;
    bool takeAbs_;
};

// Set the values of all objects in a list
void setValue(const RooArgList& args, Double_t value) {
  for (int i = 0; i < args.getSize(); ++i) {
    RooRealVar* var = dynamic_cast<RooRealVar*>(args.at(i));
    assert(var);
    var->setVal(value);
  }
}

void setError(const RooArgList& args, Double_t error) {
  for (int i = 0; i < args.getSize(); ++i) {
    RooRealVar* var = dynamic_cast<RooRealVar*>(args.at(i));
    assert(var);
    var->setError(error);
  }
}

void setConstant(const RooArgList& args, bool isConstant) {
  for (int i = 0; i < args.getSize(); ++i) {
    RooRealVar* var = dynamic_cast<RooRealVar*>(args.at(i));
    assert(var);
    var->setConstant(isConstant);
  }
}

void setRange(const RooArgList& args, double low, double high) {
  for (int i = 0; i < args.getSize(); ++i) {
    RooRealVar* var = dynamic_cast<RooRealVar*>(args.at(i));
    assert(var);
    var->setRange(low, high);
  }
}

RooArgList buildVariableCollection(size_t nVars, const std::string& namebase,
    const std::string& title, double value, double min, double max) {
  std::cout << "Creating " << nVars << " variables.  Name: " << namebase
    << std::endl;
  RooArgList output;
  for (size_t i = 0; i < nVars; ++i) {
    std::stringstream name;
    name << namebase << "_" << i;
    RooRealVar* toAdd = new RooRealVar(
        name.str().c_str(), title.c_str(), value, min, max);
    output.addOwned(*toAdd);
  }
  return output;
}

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

    double scaledMinMass = (minMass - 60.0)/100.0;
    double scaledMaxMass = (maxMass - 60.0)/100.0;

    std::stringstream massSliceCut;
    massSliceCut << scaledMinMass
      << " < scaledMass && scaledMass < " << scaledMaxMass;

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
      fittedModel.plotOn(frame.get(), ProjWData(conditionalVars, *phiSlice),
          Components("*gauss*"), LineStyle(kDotted)) ;
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
  RooFormulaVar scaledMassFormula("scaledMass", "(M-60)/100",
      "(resonanceMass-60.0)/100.0", RooArgSet(resonanceMass));

  RooRealVar* scaledLeg1Pt = static_cast<RooRealVar*>(
      data.addColumn(scaledLeg1PtFunc));
  scaledLeg1Pt->setRange(0, 15);
  scaledLeg1Pt->setBins(150);
  RooRealVar* dPhi = static_cast<RooRealVar*>(
      data.addColumn(dPhiFormula));
  dPhi->setRange(0, 1.8);
  RooRealVar* scaledMass = static_cast<RooRealVar*>(
      data.addColumn(scaledMassFormula));
  scaledMass->setRange(0, 6);

  std::cout << "Initial data set has " << data.numEntries() << std::endl;

  // Apply analysys selections and reduce the total # of variables
  RooAbsData* selectedData = data.reduce(
      RooArgSet(*scaledLeg1Pt, *dPhi, resonanceMass, *scaledMass),
      "leg1VisPt > 15 && leg2VisPt > 20"
      "&& abs(leg1VisEta) < 2.1 && abs(leg2VisEta) < 2.3");

  std::cout << "After analysis selection data has "
    << selectedData->numEntries() << "entries" <<  std::endl;

  double reducedMassMin = 245;
  double reducedMassMax = 255;

  std::stringstream reducedMassCut;
  reducedMassCut << reducedMassMax << " > resonanceMass && resonanceMass > "
    << reducedMassMin;

  RooAbsData* reducedDataSet = selectedData->reduce(
      reducedMassCut.str().c_str());

  RooAbsData* backToBackData = selectedData->reduce("dPhi < 0.03");

  std::cout << "Making reduced data set with "
    << reducedDataSet->numEntries() << " entries" <<  std::endl;

  std::cout << "Making back-to-back data set with "
    << backToBackData->numEntries() << " entries" <<  std::endl;

  RooRealVar zero("zero", "zero", 0, -1, 1);
  zero.setConstant(true);

  std::cout << "Building variable collections" << std::endl;

  RooArgList gammaScaleVars = buildVariableCollection(
      9, "GammaScale", "#Gamma scale", 0, -100, 100);

  RooArgList gammaShapeVars = buildVariableCollection(
      9, "GammaShape", "#Gamma shape", 0, -10, 10);

  RooArgList gaussMeanVars = buildVariableCollection(
      9, "GaussMean", "#Gauss mean", 0, -10, 10);

  RooArgList gaussSigmaVars = buildVariableCollection(
      9, "GaussSigma", "#Gauss sigma", 0, -10, 10);

  RooArgList gaussSigmaLVars = buildVariableCollection(
      9, "GaussSigmaL", "#Gauss sigma left", 0, -10, 10);

  RooArgList mixVars = buildVariableCollection(
      9, "Mix", "#Gauss sigma", 0, -10, 10);

  RooDoublePolyVar gammaScale("GammaScale", "#Gamma scale",
      *dPhi, *scaledMass, gammaScaleVars, 2, 2, true);

  RooDoublePolyVar gammaShape("GammaShape", "#Gamma shape",
      *dPhi, *scaledMass, gammaShapeVars, 2, 2, true);

  RooDoublePolyVar gaussMean("GaussMean", "Gaussian mean",
      *dPhi, *scaledMass, gaussMeanVars, 2, 2, true);

  RooDoublePolyVar gaussSigma("GaussSigma", "Gaussian sigma",
      *dPhi, *scaledMass, gaussSigmaVars, 2, 2, true);

  RooDoublePolyVar gaussSigmaL("GaussSigmaL", "Gaussian left sigma",
      *dPhi, *scaledMass, gaussSigmaLVars, 2, 2, true);

  RooDoublePolyVar mixRaw("MixRaw", "Mixing fraction raw",
      *dPhi, *scaledMass, mixVars, 2, 2, true);

  RooFormulaVar mix("Mix", "(1.0/3.14159)*atan(@0)+0.5", mixRaw);

  RooGamma gamma("gamma", "gamma",
      *scaledLeg1Pt, gammaScale, gammaShape, zero);

  RooGaussian gauss("gauss", "gauss",
      *scaledLeg1Pt, gaussMean, gaussSigma);

  RooBifurGauss bifurGauss("bifurGauss", "bifurGauss",
      *scaledLeg1Pt, gaussMean, gaussSigmaL, gaussSigma);

  // Build model
  RooAddPdf model("sum", "sum", gamma, bifurGauss, mix);

  model.Print("v");

  // Give good initial guess for the constant terms
  gammaScale.getTerm(0) = 14;
  gammaShape.getTerm(0) = 0.1;
  gaussMean.getTerm(0) = 1.0;
  gaussSigma.getTerm(0) = 0.1;
  gaussSigmaL.getTerm(0) = 0.1;
  mixRaw.getTerm(0) = 0.06;
  mixRaw.getTerm(1) = -1.3;

  // Taken from a previous succesfull fit.
  gammaScale.getTerm(0) = 1.40468e+01;
  gammaScale.getTerm(3) = -1.10305e+01;
  gammaScale.getTerm(6) = 3.36927e+00;
  gammaShape.getTerm(0) = 6.17678e-02;
  gammaShape.getTerm(3) = 5.82923e-02;
  gammaShape.getTerm(6) = 4.82187e-02;
  gaussMean.getTerm(0) = 9.56549e-01;
  gaussMean.getTerm(3) = -1.12297e-01;
  gaussMean.getTerm(6) = 8.45247e-02;
  gaussSigma.getTerm(0) = 7.46661e-02;
  gaussSigma.getTerm(3) = 1.92568e-01;
  gaussSigma.getTerm(6) = -3.98217e-02;
  gaussSigmaL.getTerm(0) = 7.46661e-02;
  gaussSigmaL.getTerm(3) = 1.92568e-01;
  gaussSigmaL.getTerm(6) = -3.98217e-02;
  mixRaw.getTerm(0) = 4.11155e-01;
  mixRaw.getTerm(3) = -2.74745e-01;
  mixRaw.getTerm(6) = -1.02066e+00;

  // Set the errors on all the terms to be sane, so MINUIT doesn't go crazy
  // when trying to construct the initial error matrix.
  setError(gammaScaleVars, 2e-1);
  setError(gammaShapeVars, 1e-2);
  setError(gaussMeanVars, 5e-2);
  setError(gaussSigmaVars, 5e-2);
  setError(gaussSigmaLVars, 5e-2);
  setError(mixVars, 0.3);

  // Set all the terms proportional to mass to be constant.  In the first
  // iteration we only fit the terms proportional to DPhi.
  setConstant(gammaScale.termsProportionalToY(), true);
  setConstant(gammaShape.termsProportionalToY(), true);
  setConstant(gaussMean.termsProportionalToY(), true);
  setConstant(gaussSigma.termsProportionalToY(), true);
  setConstant(gaussSigmaL.termsProportionalToY(), true);
  setConstant(mixRaw.termsProportionalToY(), true);

  std::cout << "Fitting reduced model..." << std::endl;
  RooArgSet conditionalVariables(*dPhi, *scaledMass);
  std::cout << "Begining fit" << std::endl;

  RooDataHist binnedReducedData("binnedReducedData", "binnedReducedData",
      //RooArgSet(*dPhi, *scaledLeg1Pt, resonanceMass, *scaledMass),
      RooArgSet(*dPhi, *scaledLeg1Pt, *scaledMass),
      *reducedDataSet);

  std::auto_ptr<RooPlot> scaledMassFrame(scaledMass->frame(Range(0., 10.)));
  binnedReducedData.plotOn(scaledMassFrame.get());
  scaledMassFrame->Draw();
  canvas.SaveAs("scaledMassReduced.png");

  model.fitTo(binnedReducedData, ConditionalObservables(conditionalVariables),
      NumCPU(8)
      //,Verbose(true)
      );

  std::stringstream reducedMassFitPrefix;
  reducedMassFitPrefix << "reduced_fit_mass_" << reducedMassMin
    << "_" << reducedMassMax;

  plotSlices(*scaledLeg1Pt, model, binnedReducedData, conditionalVariables,
      &canvas,reducedMassFitPrefix.str(), 1, reducedMassMin, reducedMassMax, 4, 0, 3.14/2);

  // Now do a fit to get the terms proportional to the mass.

  /*
  RooDataHist backToBackBinnedData("backToBackBinnedData",
      "backToBackBinnedData",
      RooArgSet(*dPhi, *scaledLeg1Pt, resonanceMass, *scaledMass),
      *backToBackData);
      */

  // Reduce size.
  RooDataHist backToBackBinnedData("backToBackBinnedDataForFit",
      "backToBackBinnedData",
      RooArgSet(*dPhi, *scaledLeg1Pt, *scaledMass),
      *backToBackData);

  backToBackBinnedData.plotOn(scaledMassFrame.get());
  scaledMassFrame->Draw();
  canvas.SaveAs("scaledMassBackToBack.png");

  // Now set the DPhi proportional terms constant.
  setConstant(gammaScale.termsProportionalToY(), false);
  setConstant(gammaShape.termsProportionalToY(), false);
  setConstant(gaussMean.termsProportionalToY(), false);
  setConstant(gaussSigma.termsProportionalToY(), false);
  setConstant(gaussSigmaL.termsProportionalToY(), false);
  setConstant(mixRaw.termsProportionalToY(), false);

  // And enable the terms proportional to the mass.
  setConstant(gammaScale.termsProportionalToX(), true);
  setConstant(gammaShape.termsProportionalToX(), true);
  setConstant(gaussMean.termsProportionalToX(), true);
  setConstant(gaussSigma.termsProportionalToX(), true);
  setConstant(gaussSigmaL.termsProportionalToX(), true);
  setConstant(mixRaw.termsProportionalToX(), true);

  // Don't let these blow up.
  setError(gammaScale.termsProportionalToY(), 0.01);
  setError(gaussMean.termsProportionalToY(), 0.01);

  setRange(gammaScale.termsProportionalToY(), -15, 15);
  setRange(gaussMean.termsProportionalToY(), -15, 15);

  std::auto_ptr<RooArgSet> parameters(
      model.getParameters(&backToBackBinnedData));

  std::cout << "Printing parameters." << std::endl;
  parameters->Print("v");

  std::auto_ptr<RooArgSet> observables(
      model.getObservables(&backToBackBinnedData));
  std::cout << "Printing observables." << std::endl;
  observables->Print("v");


  model.fitTo(*backToBackData,
      ConditionalObservables(conditionalVariables),
      NumCPU(8), Verbose(false), InitialHesse(false));

  plotSlices(*scaledLeg1Pt, model, backToBackBinnedData, conditionalVariables,
      &canvas, "reduced_fit_back_to_back_", 5, 100, 300, 1, 0, 3.14/2);

  RooDataHist binnedSelectedData("binnedSelectedData", "binnedSelectedData",
      RooArgSet(*dPhi, *scaledLeg1Pt, *scaledMass), *selectedData);

  plotSlices(*scaledLeg1Pt, model, binnedSelectedData, conditionalVariables,
      &canvas, "reduced_fit_all_data", 5, 100, 300, 5, 0, 3.14/2);

  setConstant(gammaScale.termsProportionalToX(), false);
  setConstant(gammaShape.termsProportionalToX(), false);
  setConstant(gaussMean.termsProportionalToX(), false);
  setConstant(gaussSigma.termsProportionalToX(), false);
  setConstant(gaussSigmaL.termsProportionalToX(), false);
  setConstant(mixRaw.termsProportionalToX(), false);

  std::cout << "Doing final fit" << std::endl;

  model.fitTo(*selectedData,
      ConditionalObservables(conditionalVariables),
      NumCPU(8), Verbose(false), InitialHesse(false));

  plotSlices(*scaledLeg1Pt, model, binnedSelectedData, conditionalVariables,
      &canvas, "final_fit_all_data", 5, 100, 300, 5, 0, 3.14/2);

  /*
  std::cout << "Creating binned dataset" << std::endl;

  RooDataHist binnedData("binnedData", "binnedData",
      RooArgSet(*dPhi, *scaledLeg1Pt, resonanceMass), *selectedData);

  plotSlices(*scaledLeg1Pt, modelDPhi, binnedData, conditionalVariables,
      &canvas, "reduced_fit", 3, 100, 300, 3, 0, 3.14/2);

  std::cout << "Now fitting full model...." << std::endl;

  model.fitTo(binnedData, ConditionalObservables(conditionalVariables),
      NumCPU(5));

  plotSlices(*scaledLeg1Pt, model, binnedData, conditionalVariables,
      &canvas, "full_fit", 5, 100, 300, 5, 0, 3.14/2);
      */

  return 0;
}
