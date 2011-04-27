#include <stdio.h>
#include <iostream>
#include <cmath>
#include <sstream>
#include <memory>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include "TGraphErrors.h"

#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooLinkedList.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooFitResult.h"

// PDFs
#include "RooGaussian.h"
#include "RooBifurGauss.h"
#include "RooGamma.h"
#include "RooMsgService.h"

//#include "TauAnalysis/CandidateTools/interface/NSVfitPtBalancePdfs.h"

using namespace RooFit;

// Simple class to hold a RooFitResult + a description and slice point.  Takes
// ownership of the RooFitResult

class PreFitResult {
  public:
    PreFitResult(RooFitResult* result, const std::string& description,
        double slicePoint):
      result_(result),
      description_(description),
      slicePoint_(slicePoint) {}

    double slicePoint() const { return slicePoint_; }
    const std::string& description() const { return description_; }
    const RooFitResult* result() const { return result_.get(); }

  private:
    boost::shared_ptr<RooFitResult> result_;
    std::string description_;
    double slicePoint_;
};

class RooDoublePolyVar : public RooAbsReal {
  public:
    RooDoublePolyVar(const char* name, const char* title,
        RooAbsReal& x, RooAbsReal& y, RooArgList& coeffs,
        int orderXLow, int orderXHigh,
        int orderYLow, int orderYHigh, bool takeAbs):
      RooAbsReal(name, title),
      x_("x", "Dependent x", this, x),
      y_("y", "Dependent y", this, y),
      coeffs_("coeffs", "List of coefficients", this),
      orderXLow_(orderXLow), orderYLow_(orderYLow),
      orderXHigh_(orderXHigh), orderYHigh_(orderYHigh),
      takeAbs_(takeAbs) {
        std::cout << "RooDoublePolyVar::ctor"
          << " order X low/high = " << orderXLow_ << "/" << orderXHigh_
          << " order Y low/high = " << orderYLow_ << "/" << orderYHigh_
          << " Ncoeffecients: " << coeffs.getSize() << std::endl;

        // Some coefficients may be unused
        size_t nCoeffs = (orderXHigh_-orderXLow_+1)*(orderYHigh_-orderYLow_+1);

        assert(nCoeffs <= (size_t)coeffs.getSize());
        for (size_t i = 0; i < nCoeffs; i++) {
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
      orderXLow_(other.orderXLow_),
      orderYLow_(other.orderYLow_),
      orderXHigh_(other.orderXHigh_),
      orderYHigh_(other.orderYHigh_),
      takeAbs_(other.takeAbs_) {}

    TObject* clone(const char* newName) const {
      return new RooDoublePolyVar(*this, newName);
    }

    Double_t evaluate() const {
      Double_t sum = 0;
      Double_t x = x_;
      Double_t y = y_;
      size_t coeffIndex = 0;
      for (int iX = orderXLow_; iX < (orderXHigh_+1); ++iX) {
        for (int iY = orderYLow_; iY < (orderYHigh_+1); ++iY) {
          const RooRealVar* coeff = dynamic_cast<const RooRealVar*>(
              coeffs_.at(coeffIndex));
          assert(coeff);
          sum += coeff->getVal()*std::pow(x, iX)*std::pow(y, iY);
          coeffIndex++;
        }
      }
      return takeAbs_ ? std::abs(sum) : sum;
    }

    RooArgList termsProportionalToX(int order) const {
      return termsProportionalToXImpl(order, false);
    }
    RooArgList termsProportionalToY(int order) const {
      return termsProportionalToYImpl(order, false);
    }

    // Get ALL terms proportional to X
    RooArgList termsProportionalToX() const {
      return termsProportionalToXImpl(-99, true);
    }
    // Get ALL terms proportional to Y
    RooArgList termsProportionalToY() const {
      return termsProportionalToYImpl(-99, true);
    }

    RooRealVar& termProportionalTo(int xOrder, int yOrder) const {
      RooArgList output;
      size_t coeffIndex = 0;
      for (int iX = orderXLow_; iX < (orderXHigh_+1); ++iX) {
        for (int iY = orderYLow_; iY < (orderYHigh_+1); ++iY) {
          RooRealVar* coeff = dynamic_cast<RooRealVar*>(
              coeffs_.at(coeffIndex));
          assert(coeff);
          coeffIndex++;
          if (iX == xOrder && iY == yOrder) {
            return *coeff;
          }
        }
      }
      assert(false);
      return *dynamic_cast<RooRealVar*>(coeffs_.at(0));
    }

    RooRealVar& getTerm(size_t i) const {
      RooRealVar* constant = dynamic_cast<RooRealVar*>(coeffs_.at(i));
      assert(constant);
      return *constant;
    }

  private:
    RooArgList termsProportionalToXImpl(int order, bool takeAll) const {
      RooArgList output;
      size_t coeffIndex = 0;
      for (int iX = orderXLow_; iX < (orderXHigh_+1); ++iX) {
        for (int iY = orderYLow_; iY < (orderYHigh_+1); ++iY) {
          RooRealVar* coeff = dynamic_cast<RooRealVar*>(
              coeffs_.at(coeffIndex));
          assert(coeff);
          coeffIndex++;
          if (takeAll && iX != 0) {
            output.add(*coeff);
          } else if (iX == order) {
            output.add(*coeff);
          }
        }
      }
      return output;
    }

    RooArgList termsProportionalToYImpl(int order, bool takeAll) const {
      RooArgList output;
      size_t coeffIndex = 0;
      for (int iX = orderXLow_; iX < (orderXHigh_+1); ++iX) {
        for (int iY = orderYLow_; iY < (orderYHigh_+1); ++iY) {
          RooRealVar* coeff = dynamic_cast<RooRealVar*>(
              coeffs_.at(coeffIndex));
          assert(coeff);
          coeffIndex++;
          if (takeAll && iY != 0) {
            output.add(*coeff);
          } else if (iY == order) {
            output.add(*coeff);
          }
        }
      }
      return output;
    }

    RooRealProxy x_;
    RooRealProxy y_;
    RooListProxy coeffs_;
    int orderXLow_;
    int orderYLow_;
    int orderXHigh_;
    int orderYHigh_;
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
    const RooArgSet& conditionalVars,
    const RooArgSet& binningVars,
    TPad* pad, const std::string& prefix,
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

    double scaledMinMass = (minMass)/100.0;
    double scaledMaxMass = (maxMass)/100.0;

    std::stringstream massSliceCut;
    massSliceCut << scaledMinMass
      << " < scaledMass && scaledMass < " << scaledMaxMass;

    std::auto_ptr<RooAbsData> massSlice(
        data.reduce(massSliceCut.str().c_str()));
    std::stringstream massRangeStr;
    massRangeStr << "(" << minMass << " - " << maxMass << ")";
    std::cout << " In mass slice " << iMass << massRangeStr.str() <<
      "there are " << massSlice->numEntries() << " entries " << std::endl;

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

      std::cout << "Binning final slice" << std::endl;
      RooDataHist binnedSlice("binnedSlice", "binnedSlice",
          binningVars, *phiSlice);

      std::auto_ptr<RooPlot> frame(toPlot.frame(Range(0, 4.0)));
      //phiSlice->plotOn(frame.get());
      binnedSlice.plotOn(frame.get());
      fittedModel.plotOn(frame.get(), ProjWData(conditionalVars, binnedSlice));
      fittedModel.plotOn(frame.get(), ProjWData(conditionalVars, binnedSlice),
          Components("*Gauss*"), LineColor(kRed)) ;
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

void plotFitResults(const std::vector<PreFitResult>& results,
    TPad* pad, const std::string& prefix) {
  if (!results.size())
    return;
  typedef std::map<std::string, TGraphErrors> GraphMap;
  GraphMap graphs;
  // Build all of the initial graphs
  const RooFitResult* firstResult = results[0].result();
  assert(firstResult);
  for (int iPar = 0; iPar < firstResult->floatParsFinal().getSize(); ++iPar){
    std::cout << "Building graph for fit parameter: ";
    std::string name = firstResult->floatParsFinal().at(iPar)->GetName();
    std::cout << name << std::endl;
    graphs[name] = TGraphErrors(results.size());
  }

  // Loop overall the mass/phi slices
  for (size_t iRes = 0; iRes < results.size(); ++iRes) {
    std::cout << "Plotting result #" << iRes << std::endl;
    const PreFitResult& result = results[iRes];
    assert(result.result());
    const RooArgList& pars = result.result()->floatParsFinal();
    for (int iPar = 0; iPar < pars.getSize(); ++iPar){
      std::string name = pars.at(iPar)->GetName();
      const RooRealVar* var = dynamic_cast<const RooRealVar*>(pars.at(iPar));
      assert(var);
      double x = result.slicePoint();
      graphs[name].SetPoint(iRes, x, var->getVal());
      graphs[name].SetPointError(iRes, 0, var->getError());
    }
  }
  BOOST_FOREACH(GraphMap::value_type& graph, graphs) {
    std::stringstream outputName;
    outputName << prefix << "_" << graph.first << ".png";
    pad->cd();
    graph.second.SetMarkerStyle(20);
    graph.second.SetMarkerSize(2);
    graph.second.Draw("ape");
    pad->SaveAs(outputName.str().c_str());
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
  RooFormulaVar scaledMassFormula("scaledMass", "(M)/100",
      "(resonanceMass)/100.0", RooArgSet(resonanceMass));

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

  RooAbsData* backToBackData = selectedData->reduce("dPhi < 0.03");
  std::cout << "Making back-to-back data set with "
    << backToBackData->numEntries() << " entries" <<  std::endl;

  RooAbsData* backToBackSinlgeMassData = backToBackData->reduce(
      "95 < resonanceMass && resonanceMass < 125");

  std::cout << "Making back-to-back signalMass data set with "
    << backToBackSinlgeMassData->numEntries() << " entries" <<  std::endl;

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

  /*
  RooDoublePolyVar gammaScale("GammaScale", "#Gamma scale",
      *dPhi, *scaledMass, gammaScaleVars, 0, 2, -1, 0, true);

  RooDoublePolyVar gammaShape("GammaShape", "#Gamma shape",
      *dPhi, *scaledMass, gammaShapeVars, 0, 2, -1, 1, true);

  RooDoublePolyVar gaussMean("GaussMean", "Gaussian mean",
      *dPhi, *scaledMass, gaussMeanVars, 0, 2, 0, 1, true);

  RooDoublePolyVar gaussSigma("GaussSigma", "Gaussian sigma",
      *dPhi, *scaledMass, gaussSigmaVars, 0, 2, 0, 1, true);

  RooDoublePolyVar gaussSigmaL("GaussSigmaL", "Gaussian left sigma",
      *dPhi, *scaledMass, gaussSigmaLVars, 0,2, 0, 1, true);

  RooDoublePolyVar mixRaw("MixRaw", "Mixing fraction raw",
      *dPhi, *scaledMass, mixVars, 0, 2, 0, 2, false);
      */
  RooDoublePolyVar gammaScale("GammaScale", "#Gamma scale",
      *dPhi, *scaledMass, gammaScaleVars, 0, 2, -1, 1, true);

  RooDoublePolyVar gammaShape("GammaShape", "#Gamma shape",
      *dPhi, *scaledMass, gammaShapeVars, 0, 2, -1, 1, true);

  RooDoublePolyVar gaussMean("GaussMean", "Gaussian mean",
      *dPhi, *scaledMass, gaussMeanVars, 0, 2, 0, 1, true);

  RooDoublePolyVar gaussSigma("GaussSigma", "Gaussian sigma",
      *dPhi, *scaledMass, gaussSigmaVars, 0, 2, 0, 1, true);

  RooDoublePolyVar gaussSigmaL("GaussSigmaL", "Gaussian left sigma",
      *dPhi, *scaledMass, gaussSigmaLVars, 0,2, 0, 1, true);

  RooDoublePolyVar mixRaw("MixRaw", "Mixing fraction raw",
      *dPhi, *scaledMass, mixVars, 0, 2, 0, 2, false);

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
  // Figured out by eye.
  gammaScale.termProportionalTo(0,0) = 14;
  gammaScale.termProportionalTo(0,-1) = 2;

  gammaShape.termProportionalTo(0, 0) = 0.06;
  gammaShape.termProportionalTo(0, -1) = -0.01;
  gammaShape.termProportionalTo(0, 1) = 0;

  gaussMean.termProportionalTo(0, 0) = 0.98;
  gaussMean.termProportionalTo(0, 1) = 0;

  gaussSigma.termProportionalTo(0, 0) = 0.08;
  gaussSigma.termProportionalTo(0, 1) = -0.01;

  gaussSigmaL.termProportionalTo(0, 0) = 0.08;
  gaussSigmaL.termProportionalTo(0, 1) = -0.01;

  mixRaw.termProportionalTo(0, 0) = 0;
  mixRaw.termProportionalTo(0, 1) = 0;

  // Get functional dependence from a previous successful fit
  setValue(gammaScale.termsProportionalToX(1), -1.10305e+01);
  setValue(gammaScale.termsProportionalToX(2), 3.36927e+00);
  setValue(gammaShape.termsProportionalToX(1), 5.82923e-02);
  setValue(gammaShape.termsProportionalToX(2), 4.82187e-02);
  setValue(gaussMean.termsProportionalToX(1), -1.12297e-01);
  setValue(gaussMean.termsProportionalToX(2), 8.45247e-02);
  setValue(gaussSigma.termsProportionalToX(1), 1.92568e-01);
  setValue(gaussSigma.termsProportionalToX(2), -3.98217e-02);
  setValue(gaussSigmaL.termsProportionalToX(1), 1.92568e-01);
  setValue(gaussSigmaL.termsProportionalToX(2), -3.98217e-02);
  setValue(mixRaw.termsProportionalToX(1), -2.74745e-01);
  setValue(mixRaw.termsProportionalToX(2), -1.02066e+00);

  // Set the errors on all the terms to be sane, so MINUIT doesn't go crazy
  // when trying to construct the initial error matrix.
  setError(gammaScaleVars, 2e-1);
  setError(gammaShapeVars, 1e-2);
  setError(gaussMeanVars, 5e-2);
  setError(gaussSigmaVars, 5e-2);
  setError(gaussSigmaLVars, 5e-2);
  setError(mixVars, 0.3);

  // Set all the terms proportional to DPhi to be constant.  In the first
  // iteration take DPhi = 0 slice and only fit the terms proportional to mass.
  setConstant(gammaScale.termsProportionalToX(), true);
  setConstant(gammaShape.termsProportionalToX(), true);
  setConstant(gaussMean.termsProportionalToX(), true);
  setConstant(gaussSigma.termsProportionalToX(), true);
  setConstant(gaussSigmaL.termsProportionalToX(), true);
  setConstant(mixRaw.termsProportionalToX(), true);

  std::cout << "Fitting reduced model..." << std::endl;
  RooArgSet conditionalVariables(*dPhi, *scaledMass);
  RooArgSet binningVars(*dPhi, *scaledLeg1Pt, *scaledMass);
  std::cout << "Beginning pre-fit" << std::endl;

  /*
  std::vector<PreFitResult> fitResults;

  size_t nMassSlices = 20;
  double lowestMass = 95;
  double highestMass = 295;
  for (size_t iMassSlice = 0; iMassSlice < nMassSlices; ++iMassSlice) {
    double reducedMassMin = lowestMass +
      iMassSlice*(highestMass-lowestMass)/nMassSlices;
    double reducedMassMax = lowestMass +
      (iMassSlice+1)*(highestMass-lowestMass)/nMassSlices;
    double nominalMass = (reducedMassMax - reducedMassMin)*0.5 + reducedMassMin;

    std::stringstream reducedMassCut;
    reducedMassCut << reducedMassMax << " > resonanceMass && resonanceMass > "
      << reducedMassMin;

    std::auto_ptr<RooAbsData> reducedDataSet(selectedData->reduce(
          binningVars,
          reducedMassCut.str().c_str()));
    std::cout << "Making mass slice data set with cut "
      << reducedMassCut.str() << " and "
      << reducedDataSet->numEntries() << " entries" <<  std::endl;

    scaledMass->setBins(20);
    dPhi->setBins(30);
    scaledLeg1Pt->setRange(0, 5);
    scaledLeg1Pt->setBins(100);

    RooDataHist binnedReducedData("binnedReducedData", "binnedReducedData",
        binningVars,
        *reducedDataSet);

    std::cout << " The binned data has " << binnedReducedData.numEntries()
      << " bins and " << binnedReducedData.sum(true) << "entries" << std::endl;

    RooFitResult* result = model.fitTo(binnedReducedData,
        ConditionalObservables(conditionalVariables),
        NumCPU(8), Save(true));

    std::stringstream reducedMassFitPrefix;
    reducedMassFitPrefix << "reduced_fit_mass_" << reducedMassMin
      << "_" << reducedMassMax;

    fitResults.push_back(PreFitResult(result, reducedMassFitPrefix.str(),
          nominalMass));

    plotSlices(*scaledLeg1Pt, model, *reducedDataSet,
        conditionalVariables, binningVars,
        &canvas,reducedMassFitPrefix.str(), 1, reducedMassMin, reducedMassMax,
        4, 0, 3.14/2);
  }

  plotFitResults(fitResults, &canvas, "dPhi_parameters_versus_mass");
  */

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

  // Reduce size.
  RooDataHist backToBackSingleMassBinnedData(
      "backToBackSingleMassBinnedDataForFit",
      "backToBackSingleMassBinnedData",
      RooArgSet(*dPhi, *scaledLeg1Pt, *scaledMass),
      *backToBackSinlgeMassData);

  /*
  backToBackBinnedData.plotOn(scaledMassFrame.get());
  scaledMassFrame->Draw();
  canvas.SaveAs("scaledMassBackToBack.png");
  */

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

  // Turn everything off to find constant terms
  setConstant(gammaScale.termsProportionalToY(), true);
  setConstant(gammaShape.termsProportionalToY(), true);
  setConstant(gaussMean.termsProportionalToY(), true);
  setConstant(gaussSigma.termsProportionalToY(), true);
  setConstant(gaussSigmaL.termsProportionalToY(), true);
  setConstant(mixRaw.termsProportionalToY(), true);

  model.fitTo(*backToBackSinlgeMassData,
      ConditionalObservables(conditionalVariables),
      NumCPU(8), Verbose(true), InitialHesse(false));

  plotSlices(*scaledLeg1Pt, model, *backToBackSinlgeMassData,
      conditionalVariables,
      binningVars,
      &canvas, "single_mass_100_back_to_back_", 1, 90, 120, 1, 0, 3.14/2);

  setConstant(gammaScale.termsProportionalToY(), false);
  setConstant(gammaShape.termsProportionalToY(), false);
  setConstant(gaussMean.termsProportionalToY(), false);
  setConstant(gaussSigma.termsProportionalToY(), false);
  setConstant(gaussSigmaL.termsProportionalToY(), false);
  setConstant(mixRaw.termsProportionalToY(), false);

  model.fitTo(backToBackBinnedData,
      ConditionalObservables(conditionalVariables),
      NumCPU(8), Verbose(false), InitialHesse(false));

  plotSlices(*scaledLeg1Pt, model, backToBackBinnedData, conditionalVariables,
      binningVars,
      &canvas, "reduced_fit_back_to_back_", 5, 100, 300, 1, 0, 3.14/2);

  RooDataHist binnedSelectedData("binnedSelectedData", "binnedSelectedData",
      RooArgSet(*dPhi, *scaledLeg1Pt, *scaledMass), *selectedData);

  plotSlices(*scaledLeg1Pt, model, binnedSelectedData, conditionalVariables,
      binningVars,
      &canvas, "reduced_fit_all_data", 5, 100, 300, 5, 0, 3.14/2);

  // Reduce the binning
  scaledLeg1Pt->setRange(0, 5);
  scaledLeg1Pt->setBins(50);
  dPhi->setBins(25);
  scaledMass->setBins(20);

  RooDataHist binnedSelectedDataForFit("binnedSelectedDataForFit", "binnedSelectedDataForFit",
      RooArgSet(*dPhi, *scaledLeg1Pt, *scaledMass), *selectedData);

  std::cout << "The binned data for fitting has "
    << binnedSelectedDataForFit.numEntries() << " bins." << std::endl;

  setConstant(gammaScale.termsProportionalToX(), false);
  setConstant(gammaShape.termsProportionalToX(), false);
  setConstant(gaussMean.termsProportionalToX(), false);
  setConstant(gaussSigma.termsProportionalToX(), false);
  setConstant(gaussSigmaL.termsProportionalToX(), false);
  setConstant(mixRaw.termsProportionalToX(), false);

  std::cout << "Doing final fit" << std::endl;

  model.fitTo(binnedSelectedDataForFit,
      ConditionalObservables(conditionalVariables),
      NumCPU(8), Verbose(false), InitialHesse(false));

  plotSlices(*scaledLeg1Pt, model, binnedSelectedData, conditionalVariables,
      binningVars,
      &canvas, "final_fit_all_data", 5, 100, 300, 5, 0, 3.14/2);

  /*
  std::cout << "Creating binned dataset" << std::endl;

  RooDataHist binnedData("binnedData", "binnedData",
      RooArgSet(*dPhi, *scaledLeg1Pt, resonanceMass), *selectedData);

  plotSlices(*scaledLeg1Pt, modelDPhi, binnedData, conditionalVariables,
      binningVars,
      &canvas, "reduced_fit", 3, 100, 300, 3, 0, 3.14/2);

  std::cout << "Now fitting full model...." << std::endl;

  model.fitTo(binnedData, ConditionalObservables(conditionalVariables),
      NumCPU(5));

  plotSlices(*scaledLeg1Pt, model, binnedData, conditionalVariables,
      binningVars,
      &canvas, "full_fit", 5, 100, 300, 5, 0, 3.14/2);
      */

  return 0;
}
