#include <stdio.h>
#include <iostream>
#include <cmath>
#include <sstream>
#include <memory>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include "TPaveLabel.h"
#include "TGraphErrors.h"

#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooLinkedList.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooDataHist.h"
#include "RooFitResult.h"

// PDFs
#include "RooGaussian.h"
#include "RooBifurGauss.h"
#include "RooGamma.h"
#include "RooMsgService.h"
#include "TauAnalysis/CandidateTools/interface/RooSkewNormal.h"

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
        int orderX, int orderY, bool takeAbs):
      RooAbsReal(name, title),
      x_("x", "Dependent x", this, x),
      y_("y", "Dependent y", this, y),
      coeffs_("coeffs", "List of coefficients", this),
      orderX_(orderX), orderY_(orderY),takeAbs_(takeAbs) {
        std::cout << "RooDoublePolyVar::ctor"
          << " order X: " << orderX_
          << " order Y: " << orderY_
          << " Ncoeffecients: " << coeffs.getSize() << std::endl;

        int nCoeffs = (orderX_+1)*(orderY_+1);
        assert(nCoeffs <= coeffs.getSize());

        for (int i = 0; i < nCoeffs; i++) {
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
      for (int iX = 0; iX < (orderX_+1); ++iX) {
        for (int iY = 0; iY < (orderY_+1); ++iY) {
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

    // Get all the terms that have some dependence on X AND Y
    RooArgList crossTerms() const {
      RooArgList output;
      size_t coeffIndex = 0;
      for (int iX = 0; iX < (orderX_+1); ++iX) {
        for (int iY = 0; iY < (orderY_+1); ++iY) {
          RooRealVar* coeff = dynamic_cast<RooRealVar*>(
              coeffs_.at(coeffIndex));
          assert(coeff);
          coeffIndex++;
          if (iX != 0 && iY != 0) {
            output.add(*coeff);
          }
        }
      }
      return output;
    }

    RooArgList allTerms() const {
      RooArgList output;
      for (int i = 0; i < coeffs_.getSize(); ++i) {
        RooRealVar* coeff = dynamic_cast<RooRealVar*>(coeffs_.at(i));
        assert(coeff);
        output.add(*coeff);
      }
      return output;
    }

    RooRealVar& termProportionalTo(int xOrder, int yOrder) const {
      RooArgList output;
      size_t coeffIndex = 0;
      for (int iX = 0; iX < (orderX_+1); ++iX) {
        for (int iY = 0; iY < (orderY_+1); ++iY) {
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

  private:
    RooArgList termsProportionalToXImpl(int order, bool takeAll) const {
      RooArgList output;
      size_t coeffIndex = 0;
      for (int iX = 0; iX < (orderX_+1); ++iX) {
        for (int iY = 0; iY < (orderY_+1); ++iY) {
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
      for (int iX = 0; iX < (orderX_+1); ++iX) {
        for (int iY = 0; iY < (orderY_+1); ++iY) {
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

  private:
    RooRealProxy x_;
    RooRealProxy y_;
    RooListProxy coeffs_;
    int orderX_;
    int orderY_;
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
    size_t phiSlices, double phiLow, double phiHigh,
    bool plotModel=true) {

  pad->cd();

  std::cout << "Plotting " << massSlices << " mass slices and "
    << phiSlices << " phi slices" << std::endl;

  TCanvas summaryCanvas("summary", "summary", 800, 600);

  std::stringstream summaryFileName;
  summaryFileName << prefix << "_summary.ps";

  std::stringstream summaryFileNameOpen;
  summaryFileNameOpen << summaryFileName.str() << "[";

  // Open the ps file
  summaryCanvas.Print(summaryFileNameOpen.str().c_str());

  //summaryCanvas.Divide(phiSlices, massSlices);

  // Maintain ownership of plots until summary plot is printed.
  boost::ptr_vector<RooPlot> summaryPlots;

  double phiSliceSize = (phiHigh - phiLow)/phiSlices;
  double massSliceSize = (massHigh - massLow)/massSlices;

  for (size_t iMass = 0; iMass < massSlices; ++iMass) {

    double minMass = massLow + iMass*massSliceSize;
    double maxMass = massLow + (iMass+1)*massSliceSize;

    double scaledMinMass = (minMass - 100.0)/100.0;
    double scaledMaxMass = (maxMass - 100.0)/100.0;

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

      //size_t subPadNumber = (iMass*phiSlices) + iPhi + 1;

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
      if (plotModel) {
        fittedModel.plotOn(frame.get(), ProjWData(conditionalVars, binnedSlice));
        fittedModel.plotOn(frame.get(), ProjWData(conditionalVars, binnedSlice),
            Components("*skew*"), LineColor(kRed)) ;
      }
      std::stringstream plotTitle;
      plotTitle << " Leg1Pt for mass " << massRangeStr.str()
        << " and #Delta#phi " << phiRangeStr.str();
      frame->SetTitle(plotTitle.str().c_str());

      // Add labels
      std::stringstream massSliceLabelTxt;
      massSliceLabelTxt << minMass << " < M < " << maxMass;
      TPaveLabel* massLabel = new TPaveLabel(
          0.6, 0.7, 0.85, 0.8, massSliceLabelTxt.str().c_str(), "");

      std::stringstream phiSliceLabelTxt;
      phiSliceLabelTxt << minPhi << " < #Delta #phi < " << maxPhi;
      TPaveLabel* phiLabel = new TPaveLabel(
          0.6, 0.8, 0.85, 0.9, phiSliceLabelTxt.str().c_str(), "");

      frame->addObject(massLabel); // takes ownership
      frame->addObject(phiLabel); // takes ownership

      // Make single plot
      pad->cd();
      frame->Draw();
      std::stringstream sliceFileName;
      sliceFileName << prefix << "_phi_"
        << iPhi << "_mass_" << iMass << ".png";
      pad->SaveAs(sliceFileName.str().c_str());

      // Make summary plot
      //summaryCanvas.cd(subPadNumber);
      summaryCanvas.cd();
      frame->Draw();
      summaryCanvas.Print(summaryFileName.str().c_str());
      // Keep ownership of the frame.
      //summaryPlots.push_back(frame);
    }
  }
  std::cout << "Closing summary canvas" << std::endl;
  std::stringstream summaryFileNameClose;
  summaryFileName  << summaryFileName.str() << "]";
  summaryCanvas.SaveAs(summaryFileNameClose.str().c_str());
  std::cout << "Done writing summary canvas" << std::endl;
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
  RooFormulaVar scaledMassFormula("scaledMass", "(M-100)/100",
      "(resonanceMass-100.0)/100.0", RooArgSet(resonanceMass));

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
      "leg1VisPt > 15 && leg2VisPt > 20 && "
      "abs(leg1VisEta) < 2.1 && abs(leg2VisEta) < 2.3");

  RooAbsData* selectedDataLeg1High = data.reduce(
      "leg1VisPt > 15 && leg2VisPt > 20 && "
      "leg2VisPt < leg1VisPt"
      //"leg1Pt > leg2Pt"
      );

  RooAbsData* selectedDataLeg2High = data.reduce(
      "leg1VisPt > 15 && leg2VisPt > 20 && "
      "leg2VisPt >= leg1VisPt"
      //"leg2Pt >= leg1Pt"
      );

  std::cout << "After analysis selection data has "
    << selectedData->numEntries() << "entries" <<  std::endl;

  RooAbsData* backToBackData = selectedData->reduce("dPhi < 0.01");
  std::cout << "Making back-to-back data set with "
    << backToBackData->numEntries() << " entries" <<  std::endl;

  RooRealVar zero("zero", "zero", 0, -1, 1);
  zero.setConstant(true);

  std::cout << "Building variable collections" << std::endl;

  RooArgList gammaScaleVars = buildVariableCollection(
      16, "GammaScale", "#Gamma scale", 0, -100, 100);

  RooArgList gammaShapeVars = buildVariableCollection(
      16, "GammaShape", "#Gamma shape", 0, -10, 10);

  RooArgList gaussMeanVars = buildVariableCollection(
      16, "GaussMean", "#Gauss mean", 0, -10, 10);

  RooArgList gaussSigmaVars = buildVariableCollection(
      16, "GaussSigma", "#Gauss sigma", 0, -10, 10);

  RooArgList gaussSigmaLVars = buildVariableCollection(
      16, "GaussSigmaL", "#Gauss sigma left", 0, -10, 10);

  RooArgList gaussSigmaBonusVars = buildVariableCollection(
      16, "GaussSigmaBonus", "#Gauss sigma", 0, -10, 10);

  RooArgList gaussSigmaLBonusVars = buildVariableCollection(
      16, "GaussSigmaLBonus", "#Gauss sigma left", 0, -10, 10);

  RooArgList mixVars = buildVariableCollection(
      16, "Mix", "#Gauss sigma", 0, -10, 10);

  RooArgList mixVarsBonus = buildVariableCollection(
      16, "MixBonus", "#Gauss sigma", 0, -10, 10);

  RooDoublePolyVar gammaScale("GammaScale", "#Gamma scale",
      *dPhi, *scaledMass, gammaScaleVars, 3, 3, true);

  // Using now for skew, true -> false
  RooDoublePolyVar gammaShape("GammaShape", "#Gamma shape",
      *dPhi, *scaledMass, gammaShapeVars, 3, 3, false);

  RooDoublePolyVar gaussMean("GaussMean", "Gaussian mean",
      *dPhi, *scaledMass, gaussMeanVars, 3, 3, true);

  RooDoublePolyVar gaussSigma("GaussSigma", "Gaussian sigma",
      *dPhi, *scaledMass, gaussSigmaVars, 3, 3, true);

  RooDoublePolyVar gaussSigmaL("GaussSigmaL", "Gaussian left sigma",
      *dPhi, *scaledMass, gaussSigmaLVars, 3, 3, true);

  RooDoublePolyVar gaussSigmaBonus("GaussSigmaBonus", "Gaussian sigma",
      *dPhi, *scaledMass, gaussSigmaBonusVars, 3, 3, true);

  RooDoublePolyVar gaussSigmaLBonus("GaussSigmaLBonus", "Gaussian left sigma",
      *dPhi, *scaledMass, gaussSigmaLBonusVars, 3, 3, true);

  RooDoublePolyVar mixRawBonus("MixRaw", "Mixing fraction raw",
      *dPhi, *scaledMass, mixVarsBonus, 3, 3, true);

  RooDoublePolyVar mixRaw("MixRaw", "Mixing fraction raw",
      *dPhi, *scaledMass, mixVars, 3, 3, true);

  RooFormulaVar mix("Mix", "(1.0/3.14159)*atan(@0)+0.5", mixRaw);
  RooFormulaVar mixBonus("Mix", "(1.0/3.14159)*atan(@0)+0.5", mixRawBonus);

  RooGamma gamma("gamma", "gamma",
      *scaledLeg1Pt, gammaScale, gammaShape, zero);

  RooGaussian gauss("gauss", "gauss",
      *scaledLeg1Pt, gaussMean, gaussSigma);

  RooBifurGauss bifurGauss("bifurGauss", "bifurGauss",
      *scaledLeg1Pt, gaussMean, gaussSigmaL, gaussSigma);

  RooBifurGauss bonusGauss("bonusGauss", "bonusGauss",
      *scaledLeg1Pt, gaussMean, gaussSigmaLBonus, gaussSigmaBonus);

  RooGenericPdf skewedGauss("skewedGauss", "skewedGauss",
      "TMath::Gaus(@0, @1, @2)*(1.0 + TMath::Erf(@3*((@0 - @1)/@2)/TMath::Sqrt(2)))",
      RooArgList(*scaledLeg1Pt, gaussMean, gaussSigma, gaussSigmaL));

  RooSkewNormal skewNormalModel("skewNormal", "skewNormal",
      *scaledLeg1Pt, gaussMean, gammaScale, gammaShape);

  // Build model
  //RooAddPdf model("sum", "sum", gamma, skewedGauss, mix);
  RooAddPdf model("sum", "sum", gauss, skewNormalModel, mix);

  //RooAddPdf doubleGauss("doubleGauss", "doubleGauss",
  //    bifurGauss, bonusGauss, mixBonus);
  // Bonus gauss version
//  RooAddPdf model("sum", "sum", gamma, doubleGauss, mix);

  model.Print("v");

  gaussMean.termProportionalTo(0,0) = 1.0;
  gaussSigma.termProportionalTo(0,0) = 7.46661e-02;
  gammaScale.termProportionalTo(0,0) = 7.46661e-02;
  gammaShape.termProportionalTo(0,0) = 7.46661e-02;
  mixRaw.termProportionalTo(0,0) = 4.11155e-01;
  mixRaw.termProportionalTo(1,0) = -2.74745e-01;
  mixRaw.termProportionalTo(2,0) = -1.02066e+00;

  /*
  // Taken from a previous succesfull fit.
  gammaScale.termProportionalTo(0,0) = 1.40468e+01;
  gammaScale.termProportionalTo(1,0) = -1.10305e+01;
  gammaScale.termProportionalTo(2,0) = 3.36927e+00;
  gammaShape.termProportionalTo(0,0) = 6.17678e-02;
  gammaShape.termProportionalTo(1,0) = 5.82923e-02;
  gammaShape.termProportionalTo(2,0) = 4.82187e-02;
  gaussMean.termProportionalTo(0,0) = 9.56549e-01;
  gaussMean.termProportionalTo(1,0) = -1.12297e-01;
  gaussMean.termProportionalTo(2,0) = 8.45247e-02;
  gaussSigma.termProportionalTo(0,0) = 7.46661e-02;
  gaussSigma.termProportionalTo(1,0) = 1.92568e-01;
  gaussSigma.termProportionalTo(2,0) = -3.98217e-02;
  gaussSigmaL.termProportionalTo(0,0) = 7.46661e-02;
  gaussSigmaL.termProportionalTo(1,0) = 1.92568e-01;
  gaussSigmaL.termProportionalTo(2,0) = -3.98217e-02;
  gaussSigmaBonus.termProportionalTo(0,0) = 7.46661e-02;
  gaussSigmaBonus.termProportionalTo(1,0) = 1.92568e-01;
  gaussSigmaBonus.termProportionalTo(2,0) = -3.98217e-02;
  gaussSigmaLBonus.termProportionalTo(0,0) = 7.46661e-02;
  gaussSigmaLBonus.termProportionalTo(1,0) = 1.92568e-01;
  gaussSigmaLBonus.termProportionalTo(2,0) = -3.98217e-02;
  mixRaw.termProportionalTo(0,0) = 4.11155e-01;
  mixRaw.termProportionalTo(1,0) = -2.74745e-01;
  mixRaw.termProportionalTo(2,0) = -1.02066e+00;
  mixRawBonus.termProportionalTo(0,0) = 0;
  mixRawBonus.termProportionalTo(1,0) = 0;
  mixRawBonus.termProportionalTo(2,0) = 0;
  */

  // Set the errors on all the terms to be sane, so MINUIT doesn't go crazy
  // when trying to construct the initial error matrix.
  setError(gammaScaleVars, 2e-1);
  setError(gammaShapeVars, 1e-2);
  setError(gaussMeanVars, 5e-2);
  setError(gaussSigmaVars, 5e-2);
  setError(gaussSigmaLVars, 5e-2);
  setError(gaussSigmaBonusVars, 5e-2);
  setError(gaussSigmaLBonusVars, 5e-2);
  setError(mixVars, 0.3);

  // Set all the terms proportional to mass to be constant.  In the first
  // iteration we only fit the terms proportional to DPhi.
  setConstant(gammaScale.termsProportionalToY(), true);
  setConstant(gammaShape.termsProportionalToY(), true);
  setConstant(gaussMean.termsProportionalToY(), true);
  setConstant(gaussSigma.termsProportionalToY(), true);
  setConstant(gaussSigmaL.termsProportionalToY(), true);
  setConstant(gaussSigmaBonus.termsProportionalToY(), true);
  setConstant(gaussSigmaLBonus.termsProportionalToY(), true);
  setConstant(mixRaw.termsProportionalToY(), true);
  setConstant(mixRawBonus.termsProportionalToY(), true);

  setConstant(gammaScale.termsProportionalToX(), true);
  setConstant(gammaShape.termsProportionalToX(), true);
  setConstant(gaussMean.termsProportionalToX(), true);
  setConstant(gaussSigma.termsProportionalToX(), true);
  setConstant(gaussSigmaL.termsProportionalToX(), true);
  setConstant(gaussSigmaBonus.termsProportionalToX(), true);
  setConstant(gaussSigmaLBonus.termsProportionalToX(), true);
  setConstant(mixRaw.termsProportionalToX(), true);
  setConstant(mixRawBonus.termsProportionalToX(), true);

  RooArgSet conditionalVariables(*dPhi, *scaledMass);
  RooArgSet binningVars(*dPhi, *scaledLeg1Pt, *scaledMass);

  std::cout << "Making summary plots..." << std::endl;

  //plotSlices(*scaledLeg1Pt, model, *selectedData,
      //conditionalVariables, binningVars,
      //&canvas, "summary_plots", 5, 95, 295,
      //4, 0, 3.14/2, false);

  /*
  plotSlices(*scaledLeg1Pt, model, *selectedDataLeg2High,
      conditionalVariables, binningVars,
      &canvas, "summary_plots_leg2_high", 5, 95, 295,
      4, 0, 3.14/2, false);

  plotSlices(*scaledLeg1Pt, model, *selectedDataLeg1High,
      conditionalVariables, binningVars,
      &canvas, "summary_plots_leg1_high", 5, 95, 295,
      4, 0, 3.14/2, false);
      */

  std::cout << "Fitting reduced model..." << std::endl;

  std::cout << "Beginning pre-fit" << std::endl;
  std::cout << " Taking mass = 0" << std::endl;

  std::auto_ptr<RooAbsData> massZeroSlice(
      selectedData->reduce("resonanceMass < 105 && resonanceMass > 95"));

  gaussMean.termProportionalTo(0,0).setConstant(false);
  gaussSigma.termProportionalTo(0,0).setConstant(false);
  gammaScale.termProportionalTo(0,0).setConstant(false);
  gammaShape.termProportionalTo(0,0).setConstant(false);
  mixRaw.termProportionalTo(0,0).setConstant(false);
  mixRaw.termProportionalTo(1,0).setConstant(false);
  mixRaw.termProportionalTo(2,0).setConstant(false);

  std::cout << "Prefit data has " << massZeroSlice->numEntries()
    << " entries" << std::endl;

  //setConstant(gaussMean.allTerms(), true);

  RooFitResult* result = model.fitTo(*massZeroSlice,
        ConditionalObservables(conditionalVariables),
        NumCPU(4), Save(true), Verbose(true));

  plotSlices(*scaledLeg1Pt, model, *massZeroSlice,
      conditionalVariables, binningVars,
      &canvas, "prefit_slice", 1, 95, 105,
      4, 0, 3.14/2);

  std::vector<PreFitResult> fitResults;

  size_t nMassSlices = 10;
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

    std::auto_ptr<RooAbsData> reducedDataSet(selectedDataLeg1High->reduce(
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

  setConstant(gaussMean.allTerms(), true);
  RooFitResult* result = model.fitTo(binnedReducedData,
      ConditionalObservables(conditionalVariables),
      NumCPU(4), Save(true));

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

  /*
  backToBackBinnedData.plotOn(scaledMassFrame.get());
  scaledMassFrame->Draw();
  canvas.SaveAs("scaledMassBackToBack.png");
  */

  // Enable the terms proportional to the mass.
  setConstant(gammaScale.termsProportionalToY(), false);
  setConstant(gammaShape.termsProportionalToY(), false);
  setConstant(gaussMean.termsProportionalToY(), false);
  setConstant(gaussSigma.termsProportionalToY(), false);
  setConstant(gaussSigmaL.termsProportionalToY(), false);
  setConstant(mixRaw.termsProportionalToY(), false);

  setConstant(gaussSigmaBonus.termsProportionalToY(), false);
  setConstant(gaussSigmaLBonus.termsProportionalToY(), false);
  setConstant(mixRawBonus.termsProportionalToY(), false);

  // Now set the DPhi proportional terms constant.
  setConstant(gammaScale.termsProportionalToX(), true);
  setConstant(gammaShape.termsProportionalToX(), true);
  setConstant(gaussMean.termsProportionalToX(), true);
  setConstant(gaussSigma.termsProportionalToX(), true);
  setConstant(gaussSigmaL.termsProportionalToX(), true);
  setConstant(mixRaw.termsProportionalToX(), true);
  setConstant(gaussSigmaBonus.termsProportionalToX(), true);
  setConstant(gaussSigmaLBonus.termsProportionalToX(), true);
  setConstant(mixRawBonus.termsProportionalToX(), true);

  // Disable the constant terms
  gammaScale.termProportionalTo(0,0).setConstant(true);
  gammaShape.termProportionalTo(0,0).setConstant(true);
  gaussMean.termProportionalTo(0,0).setConstant(true);
  gaussSigma.termProportionalTo(0,0).setConstant(true);
  gaussSigmaL.termProportionalTo(0,0).setConstant(true);
  mixRaw.termProportionalTo(0,0).setConstant(true);
  gaussSigmaBonus.termProportionalTo(0,0).setConstant(true);
  gaussSigmaLBonus.termProportionalTo(0,0).setConstant(true);
  mixRawBonus.termProportionalTo(0,0).setConstant(true);

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

  setConstant(gaussMean.allTerms(), true);
  model.fitTo(*backToBackData,
      ConditionalObservables(conditionalVariables),
      NumCPU(4), Verbose(true), InitialHesse(false));

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

  std::cout << "Fitting cross terms" << std::endl;

  // Now set everything constant, then release the cross terms
  setConstant(gammaScale.termsProportionalToY(), true);
  setConstant(gammaShape.termsProportionalToY(), true);
  setConstant(gaussMean.termsProportionalToY(), true);
  setConstant(gaussSigma.termsProportionalToY(), true);
  setConstant(gaussSigmaL.termsProportionalToY(), true);
  setConstant(mixRaw.termsProportionalToY(), true);
  setConstant(gaussSigmaBonus.termsProportionalToY(), true);
  setConstant(gaussSigmaLBonus.termsProportionalToY(), true);
  setConstant(mixRawBonus.termsProportionalToY(), true);

  setConstant(gammaScale.crossTerms(), false);
  setConstant(gammaShape.crossTerms(), false);
  setConstant(gaussMean.crossTerms(), false);
  setConstant(gaussSigma.crossTerms(), false);
  setConstant(gaussSigmaL.crossTerms(), false);
  setConstant(mixRaw.crossTerms(), false);
  setConstant(gaussSigmaBonus.crossTerms(), false);
  setConstant(gaussSigmaLBonus.crossTerms(), false);
  setConstant(mixRawBonus.crossTerms(), false);


  setConstant(gaussMean.allTerms(), true);
  model.fitTo(binnedSelectedDataForFit,
      ConditionalObservables(conditionalVariables),
      NumCPU(4), Verbose(false), InitialHesse(false));


  plotSlices(*scaledLeg1Pt, model, binnedSelectedData, conditionalVariables,
      binningVars,
      &canvas, "cross_term_fit_all_data", 5, 100, 300, 5, 0, 3.14/2);

  std::cout << "Doing final fit" << std::endl;

  // Let all terms float.
  setConstant(gammaScale.allTerms(), false);
  setConstant(gammaShape.allTerms(), false);
  setConstant(gaussMean.allTerms(), false);
  setConstant(gaussSigma.allTerms(), false);
  setConstant(gaussSigmaL.allTerms(), false);
  setConstant(mixRaw.allTerms(), false);
  setConstant(gaussSigmaBonus.allTerms(), false);
  setConstant(gaussSigmaLBonus.allTerms(), false);
  setConstant(mixRawBonus.allTerms(), false);

  setConstant(gaussMean.allTerms(), true);
  model.fitTo(binnedSelectedDataForFit,
      ConditionalObservables(conditionalVariables),
      NumCPU(4), Verbose(false), InitialHesse(false));

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
