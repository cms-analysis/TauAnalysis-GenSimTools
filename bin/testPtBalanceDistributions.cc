#include <stdio.h>
#include <iostream>
#include <iomanip>
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
#include "RooExponential.h"
#include "RooProdPdf.h"
#include "RooGamma.h"
#include "RooVoigtian.h"
#include "RooLognormal.h"
#include "RooMsgService.h"
#include "TauAnalysis/CandidateTools/interface/RooSkewNormal.h"
#include "TauAnalysis/CandidateTools/interface/RooSmearedIsotropicDecayPdf.h"
#include "TauAnalysis/CandidateTools/interface/RooBifurGeneralizedExp.h"

using namespace RooFit;

class PreFitResult {
  public:
    PreFitResult(RooFitResult* result, const std::string& description,
        double massPoint, double phiPoint, int massIndex, int phiIndex):
      result_(result),
      description_(description),
      massPoint_(massPoint),
      phiPoint_(phiPoint),
      massIndex_(massIndex),
      phiIndex_(phiIndex) {}

    double mass() const { return massPoint_; }
    double phi() const { return phiPoint_; }

    double massIndex() const { return massIndex_; }
    double phiIndex() const { return phiIndex_; }

    double getPoint(const std::string& type) const {
    if (type == "mass") {
        return mass();
      } else
        return phi();
    }

    const std::string& description() const { return description_; }
    const RooFitResult* result() const { return result_.get(); }

  private:
    boost::shared_ptr<RooFitResult> result_;
    std::string description_;
    double massPoint_;
    double phiPoint_;
    int massIndex_;
    int phiIndex_;
};

void applyFitResultToModel(const RooAbsPdf& model,
    const RooAbsData& data, const RooFitResult& result) {
  std::auto_ptr<RooArgSet> parameters(model.getParameters(data));
  // Loop over the fitted parameters
  for (int iPar = 0; iPar < result.floatParsFinal().getSize(); ++iPar){
    std::cout << "Setting result for fit parameter: ";
    std::string name = result.floatParsFinal().at(iPar)->GetName();
    std::cout << name << std::endl;
    RooRealVar* var = dynamic_cast<RooRealVar*>(
        result.floatParsFinal().at(iPar));
    assert(var);
    if (parameters->find(var->GetName())) {
      std::cout << "The model has it!" << std::endl;
      RooAbsArg* modelvar = parameters->find(var->GetName());
      assert(modelvar);
      RooRealVar* modelRealVar = dynamic_cast<RooRealVar*>(modelvar);
      assert(modelRealVar);
      *modelRealVar = var->getVal();
      modelRealVar->setError(var->getError()*2.0);
    }
  }
}

void plotFitResults(const std::vector<PreFitResult>& results,
    const std::string& variable, TPad* pad, const std::string& prefix) {
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
      double x = result.getPoint(variable);
      graphs[name].SetPoint(iRes, x, var->getVal());
      graphs[name].SetPointError(iRes, 0, var->getError());
    }
  }
  BOOST_FOREACH(GraphMap::value_type& graph, graphs) {
    std::stringstream outputName;
    outputName << prefix << "_" << graph.first;
    pad->cd();
    graph.second.SetMarkerStyle(20);
    graph.second.SetMarkerSize(2);
    graph.second.Draw("ape");
    std::string graphicsFileName = outputName.str() + ".png";
    pad->SaveAs(graphicsFileName.c_str());
    std::string rootFileName = outputName.str() + ".root";
    pad->SaveAs(rootFileName.c_str());
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
  scaledLeg1Pt->setRange(0, 30);
  scaledLeg1Pt->setBins(300);
  RooRealVar* dPhi = static_cast<RooRealVar*>(
      data.addColumn(dPhiFormula));
  dPhi->setRange(0, 1.8);
  dPhi->setBins(360);
  RooRealVar* scaledMass = static_cast<RooRealVar*>(
      data.addColumn(scaledMassFormula));
  scaledMass->setRange(-5, 6);

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

  RooRealVar zero("zero", "zero", 0, -1, 1);
  zero.setConstant(true);

  std::cout << "Building models..." << std::endl;

  // Voigtian
  RooRealVar voigtMean("voightMean", "zero", 1, 0, 10);
  RooRealVar voigtSigma("voigtSigma", "voigtSigma", 0.2, 0, 10);
  RooRealVar voigtWidth("voigtWidth", "voigtWidth", 0.2, 0, 10);

  RooRealVar isoSmearedLoc("isoSmearedLoc", "isoSmearedLoc", 1, 0, 10);
  RooRealVar isoSmearedSmear("isoSmearedSmear", "isoSmearedSmear", 0.05, 0, 1);

  RooSmearedIsotropicDecayPdf isoSmearedModel("isoSmeared", "isoSmeared",
      *scaledLeg1Pt, isoSmearedLoc, isoSmearedSmear);

  RooRealVar skewNormalLoc("skewNormalLoc", "skewNormalLoc", 1, 0, 10);
  skewNormalLoc.setError(0.2);
  skewNormalLoc = 0.85;
  //skewNormalLoc.setConstant(true);

  RooRealVar skewNormalScale("skewNormalScale", "skewNormalScale", 0.1, 1e-6, 10);
  RooRealVar skewNormalSkew("skewNormalSkew", "skewNormalSkew", 0.1, -30, 30);

  RooSkewNormal skewNormalModel("skewNormal", "skewNormal",
      *scaledLeg1Pt, skewNormalLoc, skewNormalScale, skewNormalSkew);

  RooRealVar logNormalLoc("logNormalLoc", "logNormalLoc", 1, 0, 10);
  logNormalLoc.setError(0.2);
  RooRealVar logNormalScale("logNormalScale", "logNormalScale", 0.1, 0, 10);

  RooLognormal logNormalModel("logNormal", "logNormal",
      *scaledLeg1Pt, logNormalLoc, logNormalScale);

  RooRealVar skewNormalPowMix("skewNormalPowMix", "skewNormalPowMix", 0.5, 0.1, 0.9);
  skewNormalPowMix.setError(0.1);
  RooRealVar skewNormalPowScale("skewNormalPowScale", "skewNormalPowScale", 0.5, 1e-6, 10);
  RooRealVar skewNormalPowSkew("skewNormalPowSkew", "skewNormalPowSkew", 0.1, -10, 10);
  RooRealVar skewNormalPowPow("skewNormalPowPow", "skewNormalPowPow", 2, 0, 10);
  skewNormalPowPow.setConstant(true);

  RooGenericPdf skewedNormalPow("skewedNorm4", "skewedNorm4",
      "("
      "(1.0-@3)*TMath::Exp(-0.5*TMath::Power((@0-@1)/@2, 2.0)) + "
      "@3*TMath::Exp(-0.5*TMath::Power((@0-@1)/@4, 4.0))"
      ")"
      "*(1.0 + TMath::Erf(@5*((@0 - @1)/@2)/TMath::Sqrt(2)))",
      RooArgList(*scaledLeg1Pt, skewNormalLoc, skewNormalScale,
        skewNormalPowMix, skewNormalScale, skewNormalSkew));

  *scaledLeg1Pt = 1.0;

  RooRealVar gaussianLoc("gaussianLoc", "gaussianLoc", 1, 0, 10);
  RooRealVar gaussianSmear("gaussianSmear", "gaussianSmear", 0.1, 0, 1);
  RooGaussian gaussianModel("gauss","gauss", *scaledLeg1Pt, skewNormalLoc, gaussianSmear);

  RooRealVar skewNormAndGaussMix("skewNormAndGaussMix", "skewNormAndGaussMix", 0.1, 0, 0.5);
  RooAddPdf skewNormPlusGaussModel("skewNormAndGauss", "skewNormAndGauss",
      gaussianModel, skewNormalModel, skewNormAndGaussMix);

  RooRealVar negativeOne("negativeOne", "negativeOne", -1., -5., 0.);
  negativeOne.setConstant(true);

  RooRealVar positiveOne("positiveOne", "positiveOne", 1.0, 0., 5.);
  positiveOne.setConstant(true);

  RooRealVar positiveOhPoint2("positiveOhPoint2", "positiveOhPoint2", 0.2, 0., 5.);
  positiveOhPoint2.setConstant(true);

  // was 0.7
  RooRealVar skewNormalRelPointiness("skewPoint", "skewPoint", 0.5, 0.01, 0.9);

  RooFormulaVar skewNormal2Scale("skewNormal2Scale", "skewNormal2Scale",
      "@0*@1", RooArgList(skewNormalRelPointiness, skewNormalScale));

  RooRealVar skewNormalRelLoc("skewNormalRelLoc","skewNormalRelLoc",
      1.0, 0.8, 1.2);

  RooFormulaVar skewNormalLoc2("skewNormalLoc2", "skewNormalRelLoc",
      "@0*@1", RooArgList(skewNormalRelLoc, skewNormalLoc));

  RooSkewNormal skewNormal2("skewNorm2", "skewNorm2",
      *scaledLeg1Pt, skewNormalLoc, skewNormal2Scale, skewNormalSkew);

  // Voigtian
  RooVoigtian voigtModel("voigt", "voigt", *scaledLeg1Pt,
      skewNormalLoc, voigtWidth, voigtSigma);

  // was 0.9
  RooRealVar doubleSkewMix("doubleSkewMix", "doubleSkewMix", 0.1, 0, 0.99);
  RooAddPdf skewNormPlusSkewNormModel("skewNormAndSkewNorm", "skewNormAndSkewNorm",
      skewNormal2, skewNormalModel, doubleSkewMix);
      //skewNormalModel, skewNormal2, doubleSkewMix);
      //voigtModel, skewNormalModel, doubleSkewMix);

  // Actually prefer skinny skew
  RooExponential preferFatSkew("fatSkewConstraint",
      "fatSkewConstraint", doubleSkewMix, negativeOne);

  RooGaussian preferSameLoc("preferSameLoc", "preferSameLoc",
      skewNormalRelLoc, positiveOne, positiveOhPoint2);

  RooProdPdf doubleSkewConstrained("doubleSkewNormConstrained",
      "doubleSkewNormConstrained", RooArgList(
        //preferSameLoc,
        preferFatSkew,
        skewNormPlusSkewNormModel));

  RooRealVar leftExp("leftExp", "leftExp", 0.001, 0, 2);
  RooRealVar leftGauss("leftGauss", "leftGauss", 0.1, 0, 2);
  RooRealVar rightExp("rightExp", "rightExp", 0.001, 0, 2);
  RooRealVar rightGauss("rightGauss", "rightGauss", 0.1, 0, 2);
  RooRealVar pivot("pivot", "pivot", 1.0, 0.1, 6);
  pivot.setConstant(true);

  RooBifurGeneralizedExp bifurGeneral(
      "bifurGen", "bifurGen",
      *scaledLeg1Pt, pivot, leftExp, leftGauss,
      rightExp, rightGauss);

  std::cout << "value is: " << skewedNormalPow.getVal() << std::endl;

  std::vector<RooAbsPdf*> models;
  //models.push_back(&skewNormalModel);
  //models.push_back(&skewNormPlusGaussModel);
  //models.push_back(&skewNormAndGaussConstrained);
  //models.push_back(&doubleSkewConstrained);
  models.push_back(&bifurGeneral);
  //models.push_back(&skewNormPlusIsoSmearModel);
  //models.push_back(&skewedNormalPow);
  //models.push_back(&logNormalModel); // worse than skew normal
  //models.push_back(&isoSmearedModel);
  //models.push_back(&voigtModel);
  //



  std::cout << "Doing slices" << std::endl;

  size_t nPhiSlices = 18;
  double endPhi = 1.8;
  size_t nMassSlices = 45;
  double startMass = 55;
  double endMass = 505;
  nPhiSlices = 3;
  nMassSlices = 2;
  startMass = 140;
  endMass = 160;

  RooRealVar& observable = *scaledLeg1Pt;

  for (size_t iType = 1; iType < 2; ++iType) {
    RooAbsData* dataToFit = (iType == 0) ? selectedDataLeg1High : selectedDataLeg2High;
    std::string dataTypeName = (iType == 0) ? "type1" : "type2";
    std::cout << " Fitting " << dataTypeName << std::endl;
    // Keyed by model name
    std::map<std::string, std::vector<PreFitResult> > fitResults;

    BOOST_FOREACH(RooAbsPdf* model, models) {
      fitResults[model->GetName()] = std::vector<PreFitResult>();
    }

    for (size_t iPhi = 0; iPhi < nPhiSlices; ++iPhi) {
      double phiSliceSize = endPhi/nPhiSlices;
      double minPhi = phiSliceSize*iPhi;
      double maxPhi = phiSliceSize*(iPhi+1);
      double nominalPhi = minPhi + (maxPhi - minPhi)*0.5;

      /*
         skewNormalLoc = 1;
         skewNormalLoc.setError(0.2);
         skewNormalScale = 0.1;
         skewNormalSkew = 2;
         */

      std::stringstream phiCut;
      phiCut << minPhi << " < dPhi && dPhi < " << maxPhi;

      std::stringstream phiRangeStr;
      phiRangeStr << "(" << minPhi << " - " << maxPhi << ")";

      std::auto_ptr<RooAbsData> phiSliceData(dataToFit->reduce(
            RooArgSet(observable, *dPhi, resonanceMass), phiCut.str().c_str()));

      std::cout << "Initialized phi slice " << iPhi
        << " (" << minPhi << " - " << maxPhi << ") with "
        << phiSliceData->numEntries() << " entries." << std::endl;

      for (size_t iMass = 0; iMass < nMassSlices; ++iMass) {
        double massSliceSize = (endMass - startMass)*1.0/nMassSlices;
        double minMass = startMass + massSliceSize*iMass;
        double maxMass = startMass + massSliceSize*(iMass+1);
        double nominalMass = minMass + (maxMass - minMass)*0.5;

        std::stringstream massCut;
        massCut << minMass << " < resonanceMass &&  resonanceMass < "
          << maxMass;

        std::auto_ptr<RooAbsData> massSliceData(phiSliceData->reduce(
              RooArgSet(observable, resonanceMass), massCut.str().c_str()));

        std::stringstream massRangeStr;
        massRangeStr << "(" << minMass << " - " << maxMass << ")";

        std::cout << "Initialized mass slice " << iMass
          << massRangeStr.str() << " with "
          << massSliceData->numEntries() << " entries." << std::endl;

        if (!massSliceData->numEntries())
          continue;

        std::cout << "Binning histogram" << std::endl;
        RooDataHist phiMassPointData("phiMassPoint", "phiMassPoint",
            RooArgSet(observable, resonanceMass), *massSliceData);

        std::cout << "Fitting all models" << std::endl;
        std::auto_ptr<RooPlot> frame(observable.frame(Range(0, 4.)));
        phiMassPointData.plotOn(frame.get());

        int iColor = 0;
        BOOST_FOREACH(RooAbsPdf* model, models) {
          iColor += 2;
          // Find the closest previous fit by distance in iPhi-iMass space
          std::vector<PreFitResult>& previousFits = fitResults[model->GetName()];
          if (previousFits.size()) {
            double bestDistance = 100000000000;
            const PreFitResult* bestFit = NULL;
            BOOST_FOREACH(const PreFitResult& previousFit, previousFits) {
              double massDiff = 1.0*iMass - 1.0*previousFit.massIndex();
              double phiDiff = 1.0*iPhi - 1.0*previousFit.phiIndex();
              double totalDiff = massDiff*massDiff + phiDiff*phiDiff;
              if (!bestFit || totalDiff < bestDistance) {
                bestFit = &previousFit;
                bestDistance = totalDiff;
              }
            }
            assert(bestFit);
            std::cout << "Applying previous fit parameters" << std::endl;
            applyFitResultToModel(*model, *massSliceData, *(bestFit->result()));
          }

          std::cout << "Fitting model: " << model->GetName() << std::endl;
          RooFitResult* result = model->fitTo(
              //phiMassPointData,
              *massSliceData,
              NumCPU(4),
              //Verbose(true),
              Save(true));

          fitResults[model->GetName()].push_back(
              PreFitResult(result, "", nominalMass, nominalPhi,
                iMass, iPhi));

          model->plotOn(frame.get(), LineColor(kRed + iColor));
          model->plotOn(frame.get(), LineColor(kBlue + iColor), Components("skewNormal"));
        }

        std::stringstream plotTitle;
        plotTitle << " Leg1Pt for mass " << massRangeStr.str()
          << " and #Delta#phi " << phiRangeStr.str();

        canvas.SetLogy(true);
        frame->SetTitle(plotTitle.str().c_str());
        frame->Draw();
        std::stringstream sliceFileName;
        sliceFileName << std::setfill('0');
        sliceFileName << "model_test_phi_" << dataTypeName << "_"
          << std::setw(3) << iPhi << "_mass_"
          << std::setw(3) << iMass;
        std::string logFilename = sliceFileName.str() + "_log.png";
        canvas.SaveAs(logFilename.c_str());
        canvas.SetLogy(false);
        canvas.Update();
        std::string linearFileName = sliceFileName.str() + "_linear.png";
        canvas.SaveAs(linearFileName.c_str());
        std::string rootFileName = sliceFileName.str() + ".root";
        canvas.SaveAs(rootFileName.c_str());
      }
    }

    // Now plot the graphs of all the variables
    BOOST_FOREACH(RooAbsPdf* model, models) {
      std::vector<PreFitResult> results = fitResults[model->GetName()];
      // Make mass slice plots
      for (size_t iMass = 0; iMass < nMassSlices; ++iMass) {
        std::vector<PreFitResult> massSliceResults;
        BOOST_FOREACH(const PreFitResult& result, results) {
          if (result.massIndex() == iMass) {
            massSliceResults.push_back(result);
          }
        }
        std::stringstream prefix;
        prefix << std::setfill('0') ;
        prefix << model->GetName()
          << "_type_" << dataTypeName
          << "_mass_slice_" << std::setw(3) << iMass;
        plotFitResults(massSliceResults, "phi", &canvas, prefix.str());
      }

      // Make phi slice plots
      for (size_t iPhi = 0; iPhi < nPhiSlices; ++iPhi) {
        std::vector<PreFitResult> phiSliceResults;
        BOOST_FOREACH(const PreFitResult& result, results) {
          if (result.phiIndex() == iPhi) {
            phiSliceResults.push_back(result);
          }
        }
        std::stringstream prefix;
        prefix << std::setfill('0');
        prefix << model->GetName()
          << "_type_" << dataTypeName
          << "_phi_slice_" << std::setw(3) << iPhi;
        plotFitResults(phiSliceResults, "mass", &canvas, prefix.str());
      }

    }
  }

  return 0;
}
