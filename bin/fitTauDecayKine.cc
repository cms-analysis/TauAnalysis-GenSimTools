
#include "TauAnalysis/GenSimTools/bin/tauDecayKineAuxFunctions.h"
#include "TauAnalysis/FittingTools/interface/TauDecayKinePdf.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooGamma.h"
#include "RooLognormal.h"
#include "RooGenericPdf.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooCmdArg.h"
#include "RooFit.h"
#include "RooLinkedList.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooPlot.h"

#include <TString.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TROOT.h>
#include <TMath.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TArrayD.h>
#include <TLegend.h>
#include <TBenchmark.h>
#include <TRandom3.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

void pdfTimingTest(RooAbsPdf* pdf, RooRealVar* mom, RooRealVar* sep)
{
  TRandom3 rnd;

  TBenchmark benchmark;
  benchmark.Start("pdfTimingTest");

  double pdfValueSum = 0.;

  const unsigned numCalls = 100000;
  for ( unsigned iCall = 0; iCall < numCalls; ++iCall ) {
    double momValue = rnd.Uniform(15., 200.);
    double sepValue = rnd.Uniform(0., 0.50);

    mom->setVal(momValue);
    sep->setVal(sepValue);

    double pdfValue = pdf->getVal();

    pdfValueSum += pdfValue;
  }

  benchmark.Stop("pdfTimingTest");

  std::cout << "<pdfTimingTest>:" << std::endl;
  std::cout << " pdf = " << pdf->GetName() << std::endl;
  std::cout << " sum(pdfValues) = " << pdfValueSum << std::endl;
  std::cout << " time for " << numCalls << " calls: cpu time = " << benchmark.GetCpuTime("pdfTimingTest") << "s," 
	    << " real time = " << benchmark.GetRealTime("pdfTimingTest") << "s" << std::endl;
}

//
//-------------------------------------------------------------------------------
//

double compHistogramRMSgtMax(TH1* histogram)
{
  int binMax = histogram->GetMaximumBin();
  
  double xMax = histogram->GetBinCenter(binMax);

  double rms = 0.;

  double sumBinContents = 0.;

  int numBins = histogram->GetNbinsX();
  for ( int iBin = binMax + 1; iBin <= numBins; ++iBin ) {
    double binContent = histogram->GetBinContent(iBin);
    double x = histogram->GetBinCenter(iBin);

    double dx = x - xMax;
    rms += binContent*dx*dx;
     
    sumBinContents += binContent;
  }

  if ( sumBinContents != 0. ) rms = TMath::Sqrt(rms/sumBinContents);

  return rms;
}

//
//-------------------------------------------------------------------------------
//

void storePrefitResults(RooRealVar* fitParameter, double momMin, double momMax, TGraphErrors* graph)
{
  std::cout << "<storePrefitResults>:" << std::endl;

  int iPoint = graph->GetN();

  double mom = 0.5*(momMin + momMax);

  graph->SetPoint(iPoint, mom, fitParameter->getVal());

  graph->SetPointError(iPoint, 0.5*TMath::Abs(momMax - momMin), fitParameter->getError());
}

void storeFitResults(std::vector<RooRealVar*>& fitParamCoeff, std::vector<double>& fitCoeffVal, std::vector<double>& fitCoeffErr)
{
  std::cout << "<storeFitResults>:" << std::endl;
//
  unsigned numFitParameter = fitParamCoeff.size();
  std::cout << " numFitParameter = " << numFitParameter << std::endl;
  
  fitCoeffVal.resize(numFitParameter);
  fitCoeffErr.resize(numFitParameter);
 
  for ( unsigned iParameter = 0; iParameter < numFitParameter; ++iParameter ) {
    fitCoeffVal[iParameter] = fitParamCoeff[iParameter]->getVal();
    fitCoeffErr[iParameter] = fitParamCoeff[iParameter]->getError();
    std::cout << "iParameter = " << iParameter << " (" << fitParamCoeff[iParameter]->GetName() << "): " 
	      << fitCoeffVal[iParameter] << " +/- " << fitCoeffErr[iParameter] << std::endl;
  }
} 

template <typename T>
void clearCollection(std::vector<T*>& collection)
{
  for ( typename std::vector<T*>::iterator it = collection.begin();
	it != collection.end(); ++it ) {
    delete (*it);
  }
}

struct fitManager
{
  fitManager(RooRealVar* mom, RooRealVar* dR, int decayMode, const TString& label, const TArrayD& momBinning)
    : mom_(mom),
      dR_(dR),
      decayMode_(decayMode),
      label_(label),
      momBinning_(momBinning),
      prefitParamLandauMP_(momBinning.GetSize() - 1),
      prefitParamLandauWidth_(momBinning.GetSize() - 1),
      prefitParamLandauScaleFactor_(momBinning.GetSize() - 1),
      prefitLandauScale_(momBinning.GetSize() - 1),
      prefitLandau_(momBinning.GetSize() - 1), 
      prefitParamGaussianMean_(momBinning.GetSize() - 1),
      prefitParamGaussianSigma_(momBinning.GetSize() - 1),
      prefitGaussian_(momBinning.GetSize() - 1), 
      prefitParamMix_(momBinning.GetSize() - 1),
      prefitModel_(momBinning.GetSize() - 1),
      fitParamCoeffLandauMP_(5),
      fitParamLandauMP_(0),
      fitParamCoeffLandauWidth_(5),
      fitParamLandauWidth_(0),  
      fitParamCoeffLandauScaleFactor_(5),
      fitParamLandauScaleFactor_(0),  
      fitLandauScale_(0),
      fitLandau_(0),
      fitParamCoeffGaussianMean_(5),
      fitParamGaussianMean_(0),
      fitParamCoeffGaussianSigma_(5),
      fitParamGaussianSigma_(0),
      fitGaussian_(0),
      fitParamCoeffMix_(5),
      fitParamMix_(0),
      fitModel_(0)
  {
    std::cout << "<fitManager::fitManager>:" << std::endl;
    std::cout << " mom = " << mom_->GetName() << std::endl;
    std::cout << " dR = " << dR_->GetName() << std::endl;
    std::cout << " decayMode = " << getDecayMode_string(decayMode_) << std::endl;

    prefitResultLandauMP_ = new TGraphErrors();
    prefitResultLandauWidth_ = new TGraphErrors();
    prefitResultLandauScaleFactor_ = new TGraphErrors();
    prefitResultGaussianMean_ = new TGraphErrors();
    prefitResultGaussianSigma_ = new TGraphErrors();
    prefitResultMix_ = new TGraphErrors();

//--- parametrize pt/energy dependence of fitParameters by (orthogonal) Chebyshev polynomials
//   ( cf. http://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html )
    TString momParametrizationTFormula = 
      "(1./(x*x*x*x))*([0] + [1]*x + [2]*(2.0*x*x - 1.0)" 
      " + [3]*(4.0*x*x*x - 3.0*x) + [4]*(8.0*x*x*x*x - 8.0*x*x + 1.0))";
    momParametrization_forTFormula_ = TString("TMath::Abs(").Append(momParametrizationTFormula).Append(")");
    momParametrizationMix_forTFormula_ = TString("0.5*(1.0 + TMath::TanH(").Append(momParametrization_forTFormula_).Append("))");
    TString momParametrizationRooFit = 
      "(1./(@5*@5*@5*@5))*(@0 + @1*@5 + @2*(2.0*@5*@5 - 1.0)" 
      " + @3*(4.0*@5*@5*@5 - 3.0*@5) + @4*(8.0*@5*@5*@5*@5 - 8.0*@5*@5 + 1.0))";
    momParametrization_forRooFit_ = TString("TMath::Abs(").Append(momParametrizationRooFit).Append(")");
    momParametrizationMix_forRooFit_ = TString("0.5*(1.0 + TMath::TanH(").Append(momParametrization_forRooFit_).Append("))");
    
    momParamLandauMP_ = bookTF1("LandauMP", momParametrization_forTFormula_);
    momParamLandauWidth_ = bookTF1("LandauWidth", momParametrization_forTFormula_);
    momParamLandauScaleFactor_ = bookTF1("LandauScaleFactor", momParametrization_forTFormula_);
    momParamGaussianMean_ = bookTF1("Gaussian_mean", momParametrization_forTFormula_);
    momParamGaussianSigma_ = bookTF1("Gaussian_sigma", momParametrization_forTFormula_);
    momParamMix_ = bookTF1("mix", momParametrizationMix_forTFormula_);
  }

  ~fitManager()
  {
    clearCollection(prefitParamLandauMP_);
    clearCollection(prefitParamLandauWidth_);
    clearCollection(prefitParamLandauScaleFactor_);
    clearCollection(prefitLandauScale_);
    clearCollection(prefitLandau_);
    clearCollection(prefitParamGaussianMean_);
    clearCollection(prefitParamGaussianSigma_);
    clearCollection(prefitGaussian_);
    clearCollection(prefitParamMix_);
    clearCollection(prefitModel_);

    delete prefitResultLandauMP_;
    delete prefitResultLandauWidth_;
    delete prefitResultLandauScaleFactor_;
    delete prefitResultGaussianMean_;
    delete prefitResultGaussianSigma_;

    delete momParamLandauMP_;
    delete momParamLandauWidth_;
    delete momParamLandauScaleFactor_;
    delete momParamGaussianMean_;
    delete momParamGaussianSigma_;
    delete momParamMix_;
    
    clearCollection(fitParamCoeffLandauMP_);
    delete fitParamLandauMP_;
    clearCollection(fitParamCoeffLandauWidth_);
    delete fitParamLandauWidth_;
    clearCollection(fitParamCoeffLandauScaleFactor_);
    delete fitParamLandauScaleFactor_;
    delete fitLandauScale_;
    delete fitLandau_;
    clearCollection(fitParamCoeffGaussianMean_);
    delete fitParamGaussianMean_;
    clearCollection(fitParamCoeffGaussianSigma_);
    delete fitParamGaussianSigma_;
    delete fitGaussian_;
    delete fitParamMix_;
    delete fitModel_;
  }

  TF1* bookTF1(const TString& fitParameter, const TString& formula)
  {
    std::cout << "<fitManager::bookTF1>:" << std::endl;

    TString decayMode_string = getDecayMode_string(decayMode_);
    TString tf1Name = Form("%s_%s_%s_%s", fitParameter.Data(), dR_->GetName(), decayMode_string.Data(), label_.Data());
    double momMin = momBinning_[0];
    double momMax = momBinning_[momBinning_.GetSize() - 1];
    TF1* tf1 = new TF1(tf1Name.Data(), formula, momMin, momMax);
    return tf1;
  }

  void fitTF1(TF1* tf1, TGraph* graph, std::vector<double>& fitParamValues, std::vector<double>& fitParamErrors)
  {
    std::cout << "<fitManager::fitTF1>:" << std::endl;

    graph->Fit(tf1, "W0");

    unsigned numFitParameter = tf1->GetNpar();
    std::cout << " numFitParameter = " << numFitParameter << std::endl;

    fitParamValues.resize(numFitParameter);
    fitParamErrors.resize(numFitParameter);

    for ( unsigned iParameter = 0; iParameter < numFitParameter; ++iParameter ) {
      fitParamValues[iParameter] = tf1->GetParameter(iParameter);
      fitParamErrors[iParameter] = tf1->GetParError(iParameter);
    }

    TCanvas* canvas = new TCanvas("canvas", "canvas", 1, 1, 800, 600);
    canvas->SetFillColor(10);
    canvas->SetBorderSize(2);

    graph->SetMarkerColor(1);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.2);
    graph->SetLineColor(1);
    graph->SetLineWidth(1);
    graph->Draw("ap");

    tf1->SetLineWidth(2);
    tf1->SetLineColor(2);
    tf1->Draw("LSAME");

    TString outputFileName = TString("plots/").Append(tf1->GetName()).Append(".eps");

    canvas->SaveAs(outputFileName.Data());
    canvas->SaveAs(outputFileName.ReplaceAll(".eps", ".root").Data());

    delete canvas;
  }

  RooRealVar* buildPrefitParameter(const TString& name, const TString& momBinName, double startValue, double min, double max)
  {
    std::cout << "<fitManager::buildPrefitParameter>:" << std::endl;

    TString decayMode_string = getDecayMode_string(decayMode_);
    TString fitParameterName = 
      Form("%s_%s_%s_%s_%s", name.Data(), momBinName.Data(), dR_->GetName(), decayMode_string.Data(), label_.Data());
    RooRealVar* fitParameter = new RooRealVar(fitParameterName.Data(), fitParameterName.Data(), startValue, min, max);
    return fitParameter;
  }

  RooConstVar* buildConstPrefitParameter(const TString& name, const TString& momBinName, double value)
  {
    std::cout << "<buildConstPrefitParameter>:" << std::endl;

    TString decayMode_string = getDecayMode_string(decayMode_);
    TString fitParameterName = 
      Form("%s_%s_%s_%s_%s", name.Data(), momBinName.Data(), dR_->GetName(), decayMode_string.Data(), label_.Data());
    RooConstVar* fitParameter = new RooConstVar(fitParameterName.Data(), fitParameterName.Data(), value);
    return fitParameter;
  }

  void printPrefitParameter(const TString& fitParameterName, RooRealVar* fitParameter)
  {
    std::cout << fitParameterName.Data() << ": " << fitParameter->getVal() << " +/- " << fitParameter->getError() << std::endl;
  }

  void runPrefit(const TString& inputFileName, const TString& inputDirName, const TString& label, const TString& outputFileName)
  {
    std::cout << "<fitManager::runPrefit>:" << std::endl;
    std::cout << " inputFileName = " << inputFileName.Data() << std::endl;
    std::cout << " inputDirName = " << inputDirName.Data() << std::endl;
    
    TString decayMode_string = getDecayMode_string(decayMode_);

    TString inputDirName_full = Form("%s/%s/%s", inputDirName.Data(), label.Data(), decayMode_string.Data());
    std::cout << " inputDirName_full = " << inputDirName_full.Data() << std::endl;

    TFile* inputFile = new TFile(inputFileName.Data());

    TCanvas* canvas = new TCanvas("canvas", "canvas", 1, 1, 800, 600);
    canvas->SetFillColor(10);
    canvas->SetBorderSize(2);
  
    for ( int iMomBin = 0; iMomBin < (momBinning_.GetSize() - 1); ++iMomBin ) {
      double momMin = momBinning_[iMomBin];
      double momMax = momBinning_[iMomBin + 1];
    
      TString momName_string  = "Mom";
      if      ( inputDirName.Contains("DeltaR") ) momName_string  = "Pt";
      else if ( inputDirName.Contains("Angle")  ) momName_string  = "Energy";

      TString histogramName = 
	Form("%s_%s%s%2.0fto%2.0f", decayMode_string.Data(), label.Data(), momName_string.Data(), momMin, momMax);
      TH1* histogram = (TH1*)inputFile->Get(TString(inputDirName_full.Data()).Append("/").Append(histogramName.Data()));
      std::cout << " histogramName = " << histogramName.Data() << ": histogram = " << histogram << std::endl;

      TString datahistName = TString(histogramName).Append("_datahist");
      RooDataHist* datahist = new RooDataHist(datahistName.Data(), datahistName.Data(), RooArgList(*dR_), histogram);
      std::cout << " datahist = " << datahist << std::endl;

      TString momBinName = Form("%s%2.0fto%2.0f", momName_string.Data(), momMin, momMax);

      double histogramMean = histogram->GetMean();
      //if ( histogramMean > 0.15 ) histogramMean = 0.15;
      std::cout << " histogramMean = " << histogramMean << std::endl;
      double histogramMax = histogram->GetBinCenter(histogram->GetMaximumBin());
      //if ( histogramMax > 0.15 ) histogramMax = 0.15;
      std::cout << " histogramMax = " << histogramMax << std::endl;
      double histogramRMS = histogram->GetRMS();
      //if ( histogramRMS  > 0.15 ) histogramRMS  = 0.15;
      std::cout << " histogramRMS = " << histogramRMS << std::endl;
      double histogramRMSgtMax = compHistogramRMSgtMax(histogram);
      std::cout << " histogramRMSgtMax = " << histogramRMSgtMax << std::endl;
      
      prefitParamLandauMP_[iMomBin] = buildPrefitParameter("LandauMP", momBinName, 0.8*histogramMean, 0., histogramMean);
      prefitParamLandauWidth_[iMomBin] = buildPrefitParameter("LandauWidth", momBinName, histogramRMS, 1.e-3, 0.50);
      prefitParamLandauScaleFactor_[iMomBin] = buildPrefitParameter("LandauScaleFactor", momBinName, 1., 0.5, 2.0);
      TString LandauName =
	Form("Landau_%s_%s_%s_%s", decayMode_string.Data(), momBinName.Data(), dR_->GetName(), label_.Data());

      //RooConstVar* dR0 = buildConstPrefitParameter("dR0", momBinName, histogramMax);
      //RooConstVar* dRwidth = buildConstPrefitParameter("dRwidth", momBinName, histogramRMS);
      //TString dRscaledName = 
      //  Form("dRscaled_%s_%s_%s_%s", decayMode_string.Data(), momBinName.Data(), dR_->GetName(), label_.Data());
      //RooFormulaVar* dRscaled = new RooFormulaVar(dRscaledName.Data(), "(@0 - @1)/@2", RooArgList(*dR_, *dR0, *dRwidth));

      RooRealVar* prefitParamNonPeakPwrCorr1 = buildPrefitParameter("NonPeakPwrCorr1", momBinName, 0., -1.e+1, +1.e+1);
      RooRealVar* prefitParamNonPeakPwrCorr2 = buildPrefitParameter("NonPeakPwrCorr2", momBinName, 0., -1.e+2, +1.e+2);
      RooRealVar* prefitParamNonPeakBulkLocation = buildPrefitParameter("NonPeakBulkLocation", momBinName, 1.2*histogramMean, 0., 0.25);
      RooRealVar* prefitParamNonPeakBulkScale = buildPrefitParameter("NonPeakBulkScale", momBinName, histogramRMS, 1.e-6, 0.25);
      RooRealVar* prefitParamNonPeakBulkShape = buildPrefitParameter("NonPeakBulkShape", momBinName, 0.5, -1.e+2, +1.e+2);

      prefitLandau_[iMomBin] = 
	new RooGenericPdf(LandauName.Data(), LandauName.Data(),
			  //"(1.0 + @1*@0 + @2*(2.0*@0*@0 - 1.0))*TMath::LogNormal(@0, @3, @4, @5)*(1.0 - @6*TMath::TanH(@7*(@0 - @8)))",
			  "(1.e-3 + TMath::Abs(1.0 + @1*@0 + @2*(2.0*@0*@0 - 1.0)))*TMath::Exp(-TMath::Exp((@0 - @3)/@4) + @5*(@0 - @3)/@4)/(@4*TMath::Gamma(@0))",
			  RooArgList(*dR_, 
				     *prefitParamNonPeakPwrCorr1, *prefitParamNonPeakPwrCorr2, 
				     *prefitParamNonPeakBulkLocation, *prefitParamNonPeakBulkScale, *prefitParamNonPeakBulkShape));

      prefitParamGaussianMean_[iMomBin] = 
	buildPrefitParameter("GaussianMean", momBinName, histogramMax, 0.8*histogramMax, 1.2*histogramMax);
      prefitParamGaussianSigma_[iMomBin] = 
	buildPrefitParameter("GaussianSigma", momBinName, histogramRMS, 0.5*histogramRMS, 2.0*histogramRMS);
      TString GaussianName = 
	Form("Gaussian_%s_%s_%s_%s", decayMode_string.Data(), momBinName.Data(), dR_->GetName(), label_.Data());
      //prefitGaussian_[iMomBin] = 
      //  new RooGaussian(GaussianName.Data(), GaussianName.Data(), *dR_,
      //		  *prefitParamGaussianMean_[iMomBin], *prefitParamGaussianSigma_[iMomBin]);

      RooRealVar* prefitParamNonPeakBulkLocationB = buildPrefitParameter("NonPeakBulkLocationB", momBinName, 0.8*histogramMean, 0., 0.25);
      RooRealVar* prefitParamNonPeakBulkScaleB = buildPrefitParameter("NonPeakBulkScaleB", momBinName, histogramRMS, 1.e-6, 0.25);
      RooRealVar* prefitParamNonPeakBulkShapeB = buildPrefitParameter("NonPeakBulkShapeB", momBinName, 0.5, -1.e+2, +1.e+2);
      RooRealVar* prefitParamNonPeakTailScaleB = buildPrefitParameter("NonPeakTailScaleB", momBinName, 1., 0., 1.);
      RooRealVar* prefitParamNonPeakTailShapeB = buildPrefitParameter("NonPeakTailShapeB", momBinName, 1., 1.e-3, 1.e+3);
      RooRealVar* prefitParamNonPeakTailLocationB = buildPrefitParameter("NonPeakTailLocationB", momBinName, histogramMax + 2.*histogramRMS, histogramMax, 2.);

      prefitGaussian_[iMomBin] = 
	new RooGenericPdf(GaussianName.Data(), GaussianName.Data(), 
			  //"(1.0 + @1*@0 + @2*(2.0*@0*@0 - 1.0))*TMath::LogNormal(@0, @3, @4, @5)*(1.0 - @6*TMath::TanH(@7*(@0 - @8)))",
			  "TMath::Exp(-TMath::Exp(-(@0 + @1)/@2) - @3*(@0 - @1)/@2)/(@2*TMath::Gamma(@0))*(1.0 - @3*TMath::TanH(@4*(@0 - @5)))",
			  RooArgList(*dR_, 
				     *prefitParamNonPeakBulkLocationB, *prefitParamNonPeakBulkScaleB, *prefitParamNonPeakBulkShapeB,
				     *prefitParamNonPeakTailScaleB, *prefitParamNonPeakTailShapeB, *prefitParamNonPeakTailLocationB));


      prefitParamMix_[iMomBin] = buildPrefitParameter("mix", momBinName, 0.75, 0.50, 0.95);
      
      TString dRscaledName = 
        Form("dRscaled_%s_%s_%s_%s", decayMode_string.Data(), momBinName.Data(), dR_->GetName(), label_.Data());
      RooFormulaVar* dRscaled = new RooFormulaVar(dRscaledName.Data(), "@0*@1", RooArgList(*dR_, *mom_));

      RooRealVar* gmean = buildPrefitParameter("gmean", momBinName, histogramMax, 0., 12.);
      RooRealVar* g2sigma = buildPrefitParameter("g2sigma", momBinName, histogramRMS, 0., 12.);
      RooRealVar* g4sigma = buildPrefitParameter("g4sigma", momBinName, histogramRMS, 0., 12.);
      RooRealVar* C = buildPrefitParameter("C", momBinName, 0.25, 0., 1.);
      RooRealVar* mp = buildPrefitParameter("mp", momBinName, histogramMax, 0., 12.);
      RooRealVar* width = buildPrefitParameter("width", momBinName, histogramRMSgtMax, 0., 12.);
      RooRealVar* alpha = buildPrefitParameter("alpha", momBinName, 1., 1.e-2, 1.e+2);
      RooRealVar* x0 = buildPrefitParameter("x0", momBinName, histogramMax, 0.5*histogramMax, 1.5*histogramMax);
      RooRealVar* x1 = buildPrefitParameter("x1", momBinName, histogramMax + 3.*histogramRMSgtMax, histogramMax, 12.);

      TString modelName = Form("pdf_%s_%s_%s_%s", decayMode_string.Data(), momBinName.Data(), dR_->GetName(), label_.Data());
      prefitModel_[iMomBin] = 
	//	new RooAddPdf(modelName.Data(), modelName.Data(), 
	//		      *prefitLandau_[iMomBin], *prefitGaussian_[iMomBin], *prefitParamMix_[iMomBin]);
	new TauDecayKinePdf(modelName.Data(), modelName.Data(), 
			    *dRscaled, *gmean, *g2sigma, *g4sigma, *C, *mp, *width, *alpha, *x0, *x1);

      RooConstVar* gmeanConstraint_value = 
        new RooConstVar("gmeanConstraint_value", "gmeanConstraint_value", histogramMax);
      RooConstVar* gmeanConstraint_sigma =
        new RooConstVar("gmeanConstraint_sigma", "gmeanConstraint_sigma", 0.2*histogramMax);
      RooGaussian* gmeanConstraint_pdf =
        new RooGaussian("gmeanConstraint_pdf", "gmeanConstraint_pdf",
			*gmean, *gmeanConstraint_value, *gmeanConstraint_sigma);
      RooConstVar* x0Constraint_value =
        new RooConstVar("x0Constraint_value", "x0Constraint_value", histogramMax);
      RooConstVar* x0Constraint_sigma =
        new RooConstVar("x0Constraint_sigma", "x0Constraint_sigma", 0.2*histogramMax);
      RooGaussian* x0Constraint_pdf =
        new RooGaussian("x0Constraint_pdf", "x0Constraint_pdf",
			*x0, *x0Constraint_value, *x0Constraint_sigma);
    
      RooLinkedList options;
      options.Add(new RooCmdArg(RooFit::Save(true)));
      //options.Add(new RooCmdArg(RooFit::PrintLevel(-1)));
      //options.Add(new RooCmdArg(RooFit::PrintEvalErrors(-1)));
      //options.Add(new RooCmdArg(RooFit::Warnings(-1)));
      options.Add(new RooCmdArg(RooFit::ExternalConstraints(RooArgSet(*gmeanConstraint_pdf, *x0Constraint_pdf))));
      
      RooFitResult* prefitResult = prefitModel_[iMomBin]->fitTo(*datahist, options);
      std::cout << " prefit status = " << prefitResult->status() << " (converged = 0)" << std::endl;
      delete prefitResult;

      //printPrefitParameter("LandauMP", prefitParamLandauMP_[iMomBin]);
      //printPrefitParameter("LandauWidth", prefitParamLandauWidth_[iMomBin]);
      //printPrefitParameter("LandauScaleFactor", prefitParamLandauScaleFactor_[iMomBin]);
/*
      printPrefitParameter("GaussianMean", prefitParamGaussianMean_[iMomBin]);
      printPrefitParameter("GaussianSigma", prefitParamGaussianSigma_[iMomBin]);
      printPrefitParameter("NonPeakPwrCorr1", prefitParamNonPeakPwrCorr1);
      printPrefitParameter("NonPeakPwrCorr2", prefitParamNonPeakPwrCorr2);
      printPrefitParameter("NonPeakBulkLocation", prefitParamNonPeakBulkLocation);
      printPrefitParameter("NonPeakBulkScale", prefitParamNonPeakBulkScale);
      printPrefitParameter("NonPeakBulkShape", prefitParamNonPeakBulkShape);
      printPrefitParameter("NonPeakBulkLocationB", prefitParamNonPeakBulkLocationB);
      printPrefitParameter("NonPeakBulkScaleB", prefitParamNonPeakBulkScaleB);
      printPrefitParameter("NonPeakBulkShapeB", prefitParamNonPeakBulkShapeB);
      printPrefitParameter("NonPeakTailScaleB", prefitParamNonPeakTailScaleB);
      printPrefitParameter("NonPeakTailShapeB", prefitParamNonPeakTailShapeB);
      printPrefitParameter("NonPeakTailLocatioB", prefitParamNonPeakTailLocationB);
      printPrefitParameter("mix", prefitParamMix_[iMomBin]);     
 */
      printPrefitParameter("gmean", gmean);
      printPrefitParameter("g2sigma", g2sigma);
      printPrefitParameter("g4sigma", g4sigma);
      printPrefitParameter("C", C);
      printPrefitParameter("mp", mp);
      printPrefitParameter("width", width);
      printPrefitParameter("alpha", alpha);
      printPrefitParameter("x0", x0);
      printPrefitParameter("x1", x1);

      storePrefitResults(prefitParamLandauMP_[iMomBin], momMin, momMax, prefitResultLandauMP_);
      storePrefitResults(prefitParamLandauWidth_[iMomBin], momMin, momMax, prefitResultLandauWidth_);
      storePrefitResults(prefitParamLandauScaleFactor_[iMomBin], momMin, momMax, prefitResultLandauScaleFactor_);
      storePrefitResults(prefitParamGaussianMean_[iMomBin], momMin, momMax, prefitResultGaussianMean_);
      storePrefitResults(prefitParamGaussianSigma_[iMomBin], momMin, momMax, prefitResultGaussianSigma_);
      storePrefitResults(prefitParamMix_[iMomBin], momMin, momMax, prefitResultMix_);
         
      canvas->Clear();
      canvas->SetLogy();
    
      TString frameTitle = Form("%s %s: P_{T} = %2.0f..%2.0f GeV", dR_->GetName(), decayMode_string.Data(), momMin, momMax);
      RooPlot* frame = dR_->frame(RooFit::Title(frameTitle.Data()), RooFit::Bins(100));
    
      datahist->plotOn(frame);

      prefitModel_[iMomBin]->plotOn(frame);
      prefitModel_[iMomBin]->plotOn(frame, RooFit::Components(*prefitLandau_[iMomBin]), RooFit::LineStyle(kDashed));
      prefitModel_[iMomBin]->plotOn(frame, RooFit::Components(*prefitGaussian_[iMomBin]), RooFit::LineStyle(kDashDotted));

      frame->Draw();

      TLegend legend(0.70, 0.72, 0.89, 0.89, "", "brNDC");
      legend.SetBorderSize(0);
      legend.SetFillColor(0);
      TObject* graphLandau = 0;
      TObject* graphGaussian      = 0;
      TObject* graphSum           = 0;
      TObject* graphData          = 0;
      for ( int iItem = 0; iItem < frame->numItems(); ++iItem ) {
	TString itemName = frame->nameOf(iItem);
	TObject* item = frame->findObject(itemName.Data());
	if      ( itemName.Contains("smearedLandau") ) graphLandau = item;
	else if ( itemName.Contains("Gaussian")      ) graphGaussian      = item;
	else if ( itemName.Contains("pdf")           ) graphSum           = item;
	else if ( itemName.Contains("datahist")      ) graphData          = item; 
      }
      if ( graphLandau   != 0 ) legend.AddEntry(graphLandau,   "Landau",   "l");
      if ( graphGaussian != 0 ) legend.AddEntry(graphGaussian, "Gaussian", "l");
      if ( graphSum      != 0 ) legend.AddEntry(graphSum,      "Sum",      "l");
      if ( graphData     != 0 ) legend.AddEntry(graphData,     "TAUOLA",   "p");
      legend.Draw();

      TString outputFileName_i = outputFileName;
      outputFileName_i = outputFileName_i.ReplaceAll(".", TString("_").Append(momBinName).Append("."));
      
      canvas->Update();
  
      canvas->SaveAs(outputFileName_i.Data());

      delete datahist;
    }

    delete canvas;

    delete inputFile;
    
    fitTF1(momParamLandauMP_, prefitResultLandauMP_, prefitCoeffValLandauMP_, prefitCoeffErrLandauMP_);
    fitTF1(momParamLandauWidth_, prefitResultLandauWidth_, prefitCoeffValLandauWidth_, prefitCoeffErrLandauWidth_);
    fitTF1(momParamLandauScaleFactor_, prefitResultLandauScaleFactor_, prefitCoeffValLandauScaleFactor_, prefitCoeffErrLandauScaleFactor_);
    fitTF1(momParamGaussianMean_, prefitResultGaussianMean_, prefitCoeffValGaussianMean_, prefitCoeffErrGaussianMean_);
    fitTF1(momParamGaussianSigma_, prefitResultGaussianSigma_, prefitCoeffValGaussianSigma_, prefitCoeffErrGaussianSigma_);
    fitTF1(momParamMix_, prefitResultMix_, prefitCoeffValMix_, prefitCoeffErrMix_);
  }

  RooAbsReal* buildMomDependentFitParameter(const TString& name, std::vector<RooRealVar*>& fitParameterCoefficients, 
					    std::vector<double>& prefitParamCoeffValues, std::vector<double>& prefitParamCoeffErrors,
					    bool isMix = false)
  {
    std::cout << "<fitManager::buildMomDependentFitParameter>:" << std::endl;

    TString decayMode_string = getDecayMode_string(decayMode_);

    TObjArray dependents;
  
    unsigned numCoefficients = prefitParamCoeffValues.size();
    std::cout << " numCoefficients = " << numCoefficients << std::endl;

    for ( unsigned iCoeff = 0; iCoeff < numCoefficients; ++iCoeff ) {
      TString fitParameterCoeffName = 
	Form("%s_%s_%s_%s_%s_p%u", name.Data(), mom_->GetName(), dR_->GetName(), decayMode_string.Data(), label_.Data(), iCoeff);
      std::cout << "--> creating fitParameterCoeff: name = " << fitParameterCoeffName.Data() << "," 
		<< " startValue = " << prefitParamCoeffValues[iCoeff] << " +/- " << prefitParamCoeffErrors[iCoeff] << std::endl;
      RooRealVar* fitParameterCoeff = 
	new RooRealVar(fitParameterCoeffName.Data(), fitParameterCoeffName.Data(), prefitParamCoeffValues[iCoeff]);
      dependents.Add(fitParameterCoeff);
      fitParameterCoefficients[iCoeff] = fitParameterCoeff;
    }

    dependents.Add(mom_);
  
    TString fitParameterName = 
      Form("%s_%s_%s_%s_%s", name.Data(), mom_->GetName(), dR_->GetName(), decayMode_string.Data(), label_.Data());
    TString fitParameterFormula = ( isMix ) ? momParametrizationMix_forRooFit_ : momParametrization_forRooFit_;
    std::cout << "--> creating momentum-dependent fitParameter: name = " << fitParameterName.Data() << ","
	      << " formula = " << fitParameterFormula.Data() << std::endl;
    RooAbsReal* fitParameter = new RooFormulaVar(fitParameterName.Data(), fitParameterFormula.Data(), RooArgList(dependents));
    
    return fitParameter;
  }

  void runFit(RooAbsData* dataset, const TString& outputFileName)
  {
    std::cout << "<fitManager::runFit>:" << std::endl;

    TString decayMode_string = getDecayMode_string(decayMode_);
  
    RooDataHist dataset_binned("dataset_binned", "dataset_binned", RooArgSet(*mom_, *dR_), *dataset);

    fitParamLandauMP_ = 
      buildMomDependentFitParameter("LandauMP", fitParamCoeffLandauMP_,
				    prefitCoeffValLandauMP_, prefitCoeffErrLandauMP_);
    fitParamLandauWidth_ = 
      buildMomDependentFitParameter("LandauWidth", fitParamCoeffLandauWidth_,
				    prefitCoeffValLandauWidth_, prefitCoeffErrLandauWidth_);
    fitParamLandauScaleFactor_ =
      buildMomDependentFitParameter("LandauScaleFactor", fitParamCoeffLandauScaleFactor_,
				    prefitCoeffValLandauScaleFactor_, prefitCoeffErrLandauScaleFactor_);
    TString LandauScaleName =
      Form("LandauScale_%s_%s_%s_%s", decayMode_string.Data(), "AllMom", dR_->GetName(), label_.Data());
    fitLandauScale_ = 
      new RooFormulaVar(LandauScaleName.Data(), LandauScaleName.Data(), "@0/@1", RooArgSet(*dR_, *fitParamLandauScaleFactor_));
    TString LandauName =
      Form("Landau_%s_%s_%s_%s", decayMode_string.Data(), "AllMom", dR_->GetName(), label_.Data());
    fitLandau_ = 
      new RooLandau(LandauName.Data(), LandauName.Data(), *fitLandauScale_, 
		    *fitParamLandauMP_, *fitParamLandauWidth_);
    
    fitParamGaussianMean_ = 
      buildMomDependentFitParameter("GaussianMean", fitParamCoeffGaussianMean_,
				    prefitCoeffValGaussianMean_, prefitCoeffErrGaussianMean_);
    fitParamGaussianSigma_ = 
      buildMomDependentFitParameter("GaussianSigma", fitParamCoeffGaussianSigma_,
				    prefitCoeffValGaussianSigma_, prefitCoeffErrGaussianSigma_);
    TString GaussianName = 
      Form("Gaussian_%s_%s_%s_%s", decayMode_string.Data(), "AllMom", dR_->GetName(), label_.Data());
    fitGaussian_ =
      new RooGaussian(GaussianName.Data(), GaussianName.Data(), *dR_,
		      *fitParamGaussianMean_, *fitParamGaussianSigma_);
    
    fitParamMix_ = 
      buildMomDependentFitParameter("mix", fitParamCoeffMix_,
				    prefitCoeffValMix_, prefitCoeffErrMix_, true);
			
    TString modelName = Form("pdf_%s_%s_%s_%s", decayMode_string.Data(), "AllMom", dR_->GetName(), label_.Data());
    fitModel_ = 
      new RooAddPdf(modelName.Data(), modelName.Data(), 
		    *fitLandau_, *fitGaussian_, *fitParamMix_);
  
    std::cout << "--> estimate PDF timing..." << std::endl;
    pdfTimingTest(fitLandau_, mom_, dR_);
    pdfTimingTest(fitGaussian_, mom_, dR_);
    pdfTimingTest(fitModel_, mom_, dR_);
    std::cout << " done." << std::endl;
    
    RooLinkedList options;
    options.Add(new RooCmdArg(RooFit::ConditionalObservables(*mom_)));
    options.Add(new RooCmdArg(RooFit::Save(true)));
    //options.Add(new RooCmdArg(RooFit::PrintLevel(-1)));
    //options.Add(new RooCmdArg(RooFit::PrintEvalErrors(-1)));
    //options.Add(new RooCmdArg(RooFit::Warnings(-1)));

    std::cout << "--> starting fit..." << std::endl;

    //RooFitResult* fitResult = fitModel_->fitTo(*dataset, options);
    RooFitResult* fitResult = fitModel_->fitTo(dataset_binned, options); 
    std::cout << " fit status = " << fitResult->status() << " (converged = 0)" << std::endl;
    delete fitResult;    

    std::cout << " done." << std::endl;

    std::cout << "--> saving fit results..." << std::endl;
    RooWorkspace* ws = new RooWorkspace("ws", "workspace");
    ws->import(*fitModel_);
    TString wsOutputFileName = outputFileName;
    wsOutputFileName.ReplaceAll("plots/fitTauDecayKinePlots", "../../CandidateTools/data/mcTauDecayKine");
    wsOutputFileName.ReplaceAll(".eps", "_ws.root");
    ws->writeToFile(wsOutputFileName.Data());
    std::cout << " done." << std::endl;

    std::cout << "--> making control plots..." << std::endl;
    
    storeFitResults(fitParamCoeffLandauMP_, fitCoeffValLandauMP_, fitCoeffErrLandauMP_);
    storeFitResults(fitParamCoeffLandauWidth_, fitCoeffValLandauWidth_, fitCoeffErrLandauWidth_);
    storeFitResults(fitParamCoeffLandauScaleFactor_, fitCoeffValLandauScaleFactor_, fitCoeffErrLandauScaleFactor_);
    storeFitResults(fitParamCoeffGaussianMean_, fitCoeffValGaussianMean_, fitCoeffErrGaussianMean_);
    storeFitResults(fitParamCoeffGaussianSigma_, fitCoeffValGaussianSigma_, fitCoeffErrGaussianSigma_);
    storeFitResults(fitParamCoeffMix_, fitCoeffValMix_, fitCoeffErrMix_); 
    
    TCanvas* canvas = new TCanvas("canvas", "canvas", 1, 1, 600, 900);
    canvas->SetFillColor(10);
    canvas->SetBorderSize(2);

    for ( int iMomBin = 0; iMomBin < (momBinning_.GetSize() - 1); ++iMomBin ) {
      double momMin = momBinning_[iMomBin];
      double momMax = momBinning_[iMomBin + 1];

      if ( (iMomBin % 6) == 0 ) {
	canvas->Clear();
	canvas->Divide(2,3);
      }

      TVirtualPad* pad = canvas->cd((iMomBin % 6) + 1);
      pad->SetLogy();

      TString frameTitle = Form("%s %s: P_{T} = %2.0f..%2.0f GeV", dR_->GetName(), decayMode_string.Data(), momMin, momMax);
      RooPlot* frame = dR_->frame(RooFit::Title(frameTitle.Data()), RooFit::Bins(100));
      
      TString momSliceName = Form("momSlice%2.0fto%2.0f", momMin, momMax);
      std::stringstream momSliceCut;
      momSliceCut << momMin << " < " << mom_->GetName() << " && " << mom_->GetName() << " < " << momMax;
      
      RooAbsData* momSliceDataset = dataset->reduce(momSliceCut.str().data());
      std::cout << "Selected " << momSliceDataset->numEntries() << " TTree entries" 
		<< " in slice " << Form("mom = %2.0f..%2.0f GeV", momMin, momMax) << "..." << std::endl;
      RooDataHist momSliceDataset_binned("momSliceDataset_binned", "momSliceDataset_binned", RooArgSet(*mom_, *dR_), *momSliceDataset);
      momSliceDataset_binned.plotOn(frame);
     
      fitModel_->plotOn(frame, RooFit::ProjWData(RooArgSet(*mom_), momSliceDataset_binned));
      fitModel_->plotOn(frame, RooFit::Components(*fitLandau_), RooFit::LineStyle(kDashed));
      fitModel_->plotOn(frame, RooFit::Components(*fitGaussian_), RooFit::LineStyle(kDashDotted));

      frame->Draw();

      TLegend legend(0.70, 0.72, 0.89, 0.89, "", "brNDC");
      legend.SetBorderSize(0);
      legend.SetFillColor(0);
      TObject* graphLandau = 0;
      TObject* graphGaussian      = 0;
      TObject* graphSum           = 0;
      TObject* graphData          = 0;
      for ( int iItem = 0; iItem < frame->numItems(); ++iItem ) {
	TString itemName = frame->nameOf(iItem);
	TObject* item = frame->findObject(itemName.Data());
	if      ( itemName.Contains("Landau")   ) graphLandau = item;
	else if ( itemName.Contains("Gaussian") ) graphGaussian      = item;
	else if ( itemName.Contains("pdf")      ) graphSum           = item;
	else if ( itemName.Contains("datahist") ) graphData          = item; 
      }
      if ( graphLandau   != 0 ) legend.AddEntry(graphLandau,   "Landau",   "l");
      if ( graphGaussian != 0 ) legend.AddEntry(graphGaussian, "Gaussian", "l");
      if ( graphSum      != 0 ) legend.AddEntry(graphSum,      "Sum",      "l");
      if ( graphData     != 0 ) legend.AddEntry(graphData,     "TAUOLA",   "p");
      legend.Draw();

      if ( (iMomBin % 6) == 5 ) {
	canvas->Update();
	TString outputFileName_i = outputFileName;
	outputFileName_i = outputFileName_i.ReplaceAll(".", Form("_page%u.", (iMomBin / 6) + 1));
	canvas->SaveAs(outputFileName_i.Data());
      }
    }

    std::cout << " done." << std::endl;

    delete canvas;
  }

  void writeFitResults(std::ostream& stream, const TString& psetName)
  {
    std::cout << "<fitManager::writeFitResults>:" << std::endl;

    stream << psetName << " = cms.PSet(";
    writeFitParameter(stream, "LandauMP", momParametrization_forTFormula_, fitCoeffValLandauMP_);
    writeFitParameter(stream, "LandauWidth", momParametrization_forTFormula_, fitCoeffValLandauWidth_);
    writeFitParameter(stream, "LandauScaleFactor", momParametrization_forTFormula_, fitCoeffValLandauScaleFactor_);
    writeFitParameter(stream, "GaussianMean", momParametrization_forTFormula_, fitCoeffValGaussianMean_);
    writeFitParameter(stream, "GaussianSigma", momParametrization_forTFormula_, fitCoeffValGaussianSigma_);
    writeFitParameter(stream, "mix", momParametrizationMix_forTFormula_, fitCoeffValMix_);
    stream << ")" << std::endl;
  }

  void writeFitParameter(std::ostream& stream, const TString& fitParameterName, 
			 const TString& formula, const std::vector<double>& fitParamValues)
  {
    std::cout << "<fitManager::writeFitParameter>:" << std::endl;

    stream << "    " << fitParameterName << " = cms.PSet(" << std::endl;
    stream << "        formula = cms.string(" << formula.Data() << ")," << std::endl;
    stream << "        parameter = cms.vdouble(";
    for ( unsigned iParameter = 0; iParameter < fitParamValues.size(); ++iParameter ) {
      if ( iParameter > 0 ) stream << ", ";
      stream << fitParamValues[iParameter];
    }
    stream << ")" << std::endl;
    stream << "    )" << std::endl;
  } 

  RooRealVar* mom_;
  RooRealVar* dR_;
  int decayMode_;
  TString label_;
  TArrayD momBinning_;

  std::vector<RooRealVar*> prefitParamLandauMP_;
  std::vector<RooRealVar*> prefitParamLandauWidth_;
  std::vector<RooRealVar*> prefitParamLandauScaleFactor_;
  std::vector<RooFormulaVar*> prefitLandauScale_;
  std::vector<RooAbsPdf*> prefitLandau_;
  std::vector<RooRealVar*> prefitParamGaussianMean_;
  std::vector<RooRealVar*> prefitParamGaussianSigma_;
  std::vector<RooAbsPdf*> prefitGaussian_;
  std::vector<RooRealVar*> prefitParamMix_;
  std::vector<RooAbsPdf*> prefitModel_;

  TGraphErrors* prefitResultLandauMP_;
  TGraphErrors* prefitResultLandauWidth_;
  TGraphErrors* prefitResultLandauScaleFactor_;
  TGraphErrors* prefitResultGaussianMean_;
  TGraphErrors* prefitResultGaussianSigma_;
  TGraphErrors* prefitResultMix_;

  TString momParametrization_forTFormula_;
  TString momParametrizationMix_forTFormula_;
  TString momParametrization_forRooFit_;
  TString momParametrizationMix_forRooFit_;

  TF1* momParamLandauMP_;
  TF1* momParamLandauWidth_;
  TF1* momParamLandauScaleFactor_;
  TF1* momParamGaussianMean_;
  TF1* momParamGaussianSigma_;
  TF1* momParamMix_;

  std::vector<double> prefitCoeffValLandauMP_;
  std::vector<double> prefitCoeffErrLandauMP_;
  std::vector<double> prefitCoeffValLandauWidth_;
  std::vector<double> prefitCoeffErrLandauWidth_;
  std::vector<double> prefitCoeffValLandauScaleFactor_;
  std::vector<double> prefitCoeffErrLandauScaleFactor_;
  std::vector<double> prefitCoeffValGaussianMean_;
  std::vector<double> prefitCoeffErrGaussianMean_;
  std::vector<double> prefitCoeffValGaussianSigma_;
  std::vector<double> prefitCoeffErrGaussianSigma_;
  std::vector<double> prefitCoeffValMix_;
  std::vector<double> prefitCoeffErrMix_;

  
  std::vector<RooRealVar*> fitParamCoeffLandauMP_;
  RooAbsReal* fitParamLandauMP_;
  std::vector<RooRealVar*> fitParamCoeffLandauWidth_;
  RooAbsReal* fitParamLandauWidth_;
  std::vector<RooRealVar*> fitParamCoeffLandauScaleFactor_;
  RooAbsReal* fitParamLandauScaleFactor_;
  RooFormulaVar* fitLandauScale_;
  RooAbsPdf* fitLandau_;
  std::vector<RooRealVar*> fitParamCoeffGaussianMean_;
  RooAbsReal* fitParamGaussianMean_;
  std::vector<RooRealVar*> fitParamCoeffGaussianSigma_;
  RooAbsReal* fitParamGaussianSigma_;
  RooAbsPdf* fitGaussian_;
  std::vector<RooRealVar*> fitParamCoeffMix_;
  RooAbsReal*  fitParamMix_;
  RooAbsPdf* fitModel_;

  std::vector<double> fitCoeffValLandauMP_;
  std::vector<double> fitCoeffErrLandauMP_;
  std::vector<double> fitCoeffValLandauWidth_;
  std::vector<double> fitCoeffErrLandauWidth_;
  std::vector<double> fitCoeffValLandauScaleFactor_;
  std::vector<double> fitCoeffErrLandauScaleFactor_;
  std::vector<double> fitCoeffValGaussianMean_;
  std::vector<double> fitCoeffErrGaussianMean_;
  std::vector<double> fitCoeffValGaussianSigma_;
  std::vector<double> fitCoeffErrGaussianSigma_;
  std::vector<double> fitCoeffValMix_;
  std::vector<double> fitCoeffErrMix_;
};

//
//-------------------------------------------------------------------------------
//

int main(int argc, const char* argv[])
{
  if ( argc < 4 ) {
    std::cerr << "Usage: ./fitTauDecayKine inputFileNames decayMode selection" << std::endl;
    return 1;
  }

  TString inputFileNames_ntuple = argv[1];
  //inputFileNames_ntuple = "/data2/friis/PtBalanceNtupleData_v2/ptBalanceData_*.root";
  //inputFileNames_ntuple = "/data2/friis/PtBalanceNtupleData_v2/ptBalanceData_mass_200_ggAH.root";

  std::vector<unsigned> decayModesToRun;
  for ( unsigned iDecayMode = kElectron_Muon; iDecayMode <= kThreeProng1Pi0; ++iDecayMode ) {
    if ( getDecayMode_string(iDecayMode) == argv[2] ) decayModesToRun.push_back(iDecayMode);
  }

  if ( decayModesToRun.size() == 0 ) 
    throw cms::Exception("fitTauDecayKine")
      << "Invalid Configuration Parameter 'decayMode' = " << argv[2] << " !!\n";
  
  bool runAll      = ( std::string(argv[3]) == "all"      ) ? true : false;
  bool runSelected = ( std::string(argv[3]) == "selected" ) ? true : false;
  if ( !(runAll || runSelected) )
    throw cms::Exception("fitTauDecayKine")
      << "Invalid Configuration Parameter 'selection' = " << argv[3] << " !!\n";
  
  TString inputFileName_histograms = "makeTauDecayKinePlots.root";

  TArrayD momBinning = getBinningMom();
  TArrayD sepBinning = getBinningSepTimesMom();
  
  double decayModeElectron_Muon_encoded = encodeDecayMode(kElectron_Muon);
  double decayModeHadMin_encoded = encodeDecayMode(kOneProng0Pi0);
  double decayModeHadMax_encoded = encodeDecayMode(kThreeProng1Pi0);

  TString leg1DecayModeSign = ( decayModeElectron_Muon_encoded < 0. ) ? "+" : "-";
  TString leg2DecayModeSign = ( decayModeElectron_Muon_encoded < 0. ) ? "+" : "-";

  TString nanFilter;
  nanFilter.Append("!(TMath::IsNaN(leg1Pt) || TMath::IsNaN(leg1VisPt) ||");
  nanFilter.Append("  TMath::IsNaN(leg1VisInvisAngleLab) || TMath::IsNaN(leg1VisInvisDeltaRLab) ||");
  nanFilter.Append("  TMath::IsNaN(leg2Pt) || TMath::IsNaN(leg2VisPt) ||");
  nanFilter.Append("  TMath::IsNaN(leg2VisInvisAngleLab) || TMath::IsNaN(leg2VisInvisDeltaRLab))");

  TString visMomCuts;
  visMomCuts.Append(Form("((leg1VisPt > 15. && TMath::Abs(leg1VisEta) < 2.1 &&"));
  visMomCuts.Append(Form("  TMath::Abs(leg1DecayMode %s %2.1f) < 0.1) || ", 
			    leg1DecayModeSign.Data(), TMath::Abs(decayModeElectron_Muon_encoded)));
  visMomCuts.Append(Form(" (leg1VisPt > 20. && TMath::Abs(leg1VisEta) < 2.3 &&"));
  visMomCuts.Append(Form("  leg1DecayMode > %2.1f && leg1DecayMode < %2.1f))", 
			    decayModeHadMin_encoded - 0.1, decayModeHadMax_encoded + 0.1));
  visMomCuts.Append(" && ");
  visMomCuts.Append(Form("((leg2VisPt > 15. && TMath::Abs(leg2VisEta) < 2.1 &&"));
  visMomCuts.Append(Form("  TMath::Abs(leg2DecayMode %s %2.1f) < 0.1) || ", 
			    leg2DecayModeSign.Data(), TMath::Abs(decayModeElectron_Muon_encoded)));
  visMomCuts.Append(Form(" (leg2VisPt > 20. && TMath::Abs(leg2VisEta) < 2.3 &&"));
  visMomCuts.Append(Form("  leg2DecayMode > %2.1f && leg2DecayMode < %2.1f))", 
			    decayModeHadMin_encoded - 0.1, decayModeHadMax_encoded + 0.1));
  std::cout << "visMomCuts = " << visMomCuts.Data() << std::endl;
  
//--- stop RooT from keeping references to all histograms
  TH1::AddDirectory(false); 

  gROOT->SetBatch(true);

  // Load the data
  std::cout << "Loading data from " << inputFileNames_ntuple <<  std::endl;
  TChain* dataTree = new TChain("makePtBalanceNtuple/ptBalanceNtuple");
  dataTree->Add(inputFileNames_ntuple);

  /*****************************************************************************
   ********                   Define variabes and data                      ****
   *****************************************************************************/

  RooRealVar leg1Energy("leg1Energy", "leg1Energy", 0., 350.);
  RooRealVar leg1Pt("leg1Pt", "leg1Pt", 0., 250.);
  RooRealVar leg1VisPt("leg1VisPt", "leg1VisPt", 0., 200.);
  RooRealVar leg1VisEta("leg1VisEta", "leg1VisEta", -2.5, +2.5);
  RooRealVar leg1DecayMode("leg1DecayMode", "leg1DecayMode", -1.5, +11.5);
  RooRealVar leg1VisInvisAngleLab("leg1VisInvisAngleLab", "leg1VisInvisAngleLab", 0., sepBinning[sepBinning.GetSize() - 1]);
  RooRealVar leg1VisInvisDeltaRLab("leg1VisInvisDeltaRLab", "leg1VisInvisDeltaRLab", 0., sepBinning[sepBinning.GetSize() - 1]);

  leg1Energy.setBins(350);
  leg1Pt.setBins(250);
  leg1VisInvisAngleLab.setBins(120);
  leg1VisInvisDeltaRLab.setBins(120);

  RooRealVar leg2Energy("leg2Energy", "leg2Energy", 0., 350.);
  RooRealVar leg2Pt("leg2Pt", "leg2Pt", 0., 250.);
  RooRealVar leg2VisPt("leg2VisPt", "leg2VisPt", 0., 200.);
  RooRealVar leg2VisEta("leg2VisEta", "leg2VisEta", -2.5, +2.5);
  RooRealVar leg2DecayMode("leg2DecayMode", "leg2DecayMode", -1.5, +11.5);
  RooRealVar leg2VisInvisAngleLab("leg2VisInvisAngleLab", "leg2VisInvisAngleLab", 0., sepBinning[sepBinning.GetSize() - 1]);
  RooRealVar leg2VisInvisDeltaRLab("leg2VisInvisDeltaRLab", "leg2VisInvisDeltaRLab", 0., sepBinning[sepBinning.GetSize() - 1]);

  leg2Energy.setBins(350);
  leg2Pt.setBins(250);
  leg2VisInvisAngleLab.setBins(120);
  leg2VisInvisDeltaRLab.setBins(120);

  TObjArray variables;
  variables.Add(&leg1Energy);
  variables.Add(&leg1Pt);
  variables.Add(&leg1VisPt);
  variables.Add(&leg1VisEta);
  variables.Add(&leg1DecayMode);
  variables.Add(&leg1VisInvisAngleLab);
  variables.Add(&leg1VisInvisDeltaRLab);
  variables.Add(&leg2Energy);
  variables.Add(&leg2Pt);
  variables.Add(&leg2VisPt);
  variables.Add(&leg2VisEta);
  variables.Add(&leg2DecayMode);
  variables.Add(&leg2VisInvisAngleLab);
  variables.Add(&leg2VisInvisDeltaRLab);
/*
  RooAbsData* dataset_all = new RooDataSet("dataset", "datasetset", RooArgSet(variables), RooFit::Import(*dataTree));
  std::cout << "Processing " << dataset_all->numEntries() << " TTree entries..." << std::endl;
  
  RooAbsData* dataset_notNaN = dataset_all->reduce(nanFilter);
  std::cout << "--> Passing anti-NaN filter = " << dataset_notNaN->numEntries() << std::endl;
  
  RooAbsData* dataset_selected = 0;
  if ( runSelected ) {
    dataset_selected = dataset_notNaN->reduce(visMomCuts);
    std::cout << "--> Passing vis. Momentum Cuts = " << dataset_selected->numEntries() << std::endl;
  }
 */
  RooAbsData* dataset_notNaN = 0;
  RooAbsData* dataset_selected = 0;

  std::map<std::string, fitManager*> fitResults;
  
  for ( std::vector<unsigned>::const_iterator decayModeToRun = decayModesToRun.begin();
	decayModeToRun != decayModesToRun.end(); ++decayModeToRun ) {
    double decayMode_encoded = encodeDecayMode(*decayModeToRun);
  
    TString decayMode_string = getDecayMode_string(*decayModeToRun);
  
    TString leg1DecayModeSign = ( decayMode_encoded < 0. ) ? "+" : "-";
    TString leg1DecayModeSelection = Form("TMath::Abs(leg1DecayMode %s %2.1f) < 0.1", 
					  leg1DecayModeSign.Data(), TMath::Abs(decayMode_encoded));
    if ( (*decayModeToRun) == kOneProngGt0Pi0 ) 
      leg1DecayModeSelection = 
	Form("leg1DecayMode > %2.1f && leg1DecayMode < %2.1f", 
	     encodeDecayMode(kOneProng1Pi0) - 0.1, encodeDecayMode(kOneProng2Pi0) + 0.1);

    TString leg2DecayModeSign = ( decayMode_encoded < 0. ) ? "+" : "-";
    TString leg2DecayModeSelection = Form("TMath::Abs(leg2DecayMode %s %2.1f) < 0.1", 
					  leg2DecayModeSign.Data(), TMath::Abs(decayMode_encoded));
    if ( (*decayModeToRun) == kOneProngGt0Pi0 ) 
      leg2DecayModeSelection = 
	Form("leg2DecayMode > %2.1f && leg2DecayMode < %2.1f", 
	     encodeDecayMode(kOneProng1Pi0) - 0.1, encodeDecayMode(kOneProng2Pi0) + 0.1);
  
    TString selection_string = "";
    if ( runAll      ) selection_string.Append("_all");
    if ( runSelected ) selection_string.Append("_selected");
      
    /*****************************************************************************
     ******** Run fit for leg1 with no cuts on visible decay products applied ****
     *****************************************************************************/

    if ( runAll ) {
/*
      RooAbsData* dataset_leg1_all = dataset_notNaN->reduce(leg1DecayModeSelection.Data());
 */      
      RooAbsData* dataset_leg1_all = 0;

      TString fitManagerName_leg1_dR_all = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg1", "dR", "all");
      if ( fitResults.find(fitManagerName_leg1_dR_all.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg1_dR_all.Data()] = 
	  new fitManager(&leg1Pt, &leg1VisInvisDeltaRLab, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg1", "dR", "all"), momBinning);
      }
      fitManager* fitManager_leg1_dR_all = fitResults[fitManagerName_leg1_dR_all.Data()];
      TString outputFileName_leg1_dR_all = Form("plots/fitTauDecayKinePlots_%s_leg1_dR_all.eps", decayMode_string.Data());
      fitManager_leg1_dR_all->runPrefit(inputFileName_histograms, 
					"VisInvisDeltaRLabTimesPt", "all", outputFileName_leg1_dR_all);
      fitManager_leg1_dR_all->runFit(dataset_leg1_all, outputFileName_leg1_dR_all);
      
      TString fitManagerName_leg1_angle_all = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg1", "angle", "all");
      if ( fitResults.find(fitManagerName_leg1_angle_all.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg1_angle_all.Data()] = 
	  new fitManager(&leg1Energy, &leg1VisInvisAngleLab, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg1", "angle", "all"), momBinning);
      }
      fitManager* fitManager_leg1_angle_all = fitResults[fitManagerName_leg1_angle_all.Data()];
      TString outputFileName_leg1_angle_all = 
	Form("plots/fitTauDecayKinePlots_%s_leg1_angle_all.eps", decayMode_string.Data());
      fitManager_leg1_angle_all->runPrefit(inputFileName_histograms, 
					   "VisInvisAngleLabTimesEnergy", "all", outputFileName_leg1_angle_all);
      fitManager_leg1_angle_all->runFit(dataset_leg1_all, outputFileName_leg1_angle_all);
    }

    /*****************************************************************************
     ********   Run fit for leg1 with cuts on visible decay products applied  ****
     *****************************************************************************/

    if ( runSelected ) {
      RooAbsData* dataset_leg1_selected = dataset_selected->reduce(leg1DecayModeSelection.Data());
      
      TString fitManagerName_leg1_dR_selected = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg1", "dR", "selected");
      if ( fitResults.find(fitManagerName_leg1_dR_selected.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg1_dR_selected.Data()] = 
	  new fitManager(&leg1Pt, &leg1VisInvisDeltaRLab, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg1", "dR", "selected2"), momBinning);
      }
      fitManager* fitManager_leg1_dR_selected = fitResults[fitManagerName_leg1_dR_selected.Data()];
      TString outputFileName_leg1_dR_selected = Form("plots/fitTauDecayKinePlots_%s_leg1_dR_selected.eps", decayMode_string.Data());
      fitManager_leg1_dR_selected->runPrefit(inputFileName_histograms, 
					   "VisInvisDeltaRLabTimesPt", "selected2", outputFileName_leg1_dR_selected);
      fitManager_leg1_dR_selected->runFit(dataset_leg1_selected, outputFileName_leg1_dR_selected);
      
      TString fitManagerName_leg1_angle_selected = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg1", "angle", "selected");
      if ( fitResults.find(fitManagerName_leg1_angle_selected.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg1_angle_selected.Data()] = 
	  new fitManager(&leg1Energy, &leg1VisInvisAngleLab, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg1", "angle", "selected2"), momBinning);
      }
      fitManager* fitManager_leg1_angle_selected = fitResults[fitManagerName_leg1_angle_selected.Data()];
      TString outputFileName_leg1_angle_selected = 
	Form("plots/fitTauDecayKinePlots_%s_leg1_angle_selected.eps", decayMode_string.Data());
      fitManager_leg1_angle_selected->runPrefit(inputFileName_histograms, 
						"VisInvisAngleLabTimesEnergy", "selected2", outputFileName_leg1_angle_selected);
      fitManager_leg1_angle_selected->runFit(dataset_leg1_selected, outputFileName_leg1_angle_selected);
    }

    /*****************************************************************************
     ******** Run fit for leg2 with no cuts on visible decay products applied ****
     *****************************************************************************/

    if ( runAll ) {
      RooAbsData* dataset_leg2_all = dataset_notNaN->reduce(leg2DecayModeSelection.Data());
      
      TString fitManagerName_leg2_dR_all = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg2", "dR", "all");
      if ( fitResults.find(fitManagerName_leg2_dR_all.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg2_dR_all.Data()] = 
	  new fitManager(&leg2Pt, &leg2VisInvisDeltaRLab, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg2", "dR", "all"), momBinning);
      }
      fitManager* fitManager_leg2_dR_all = fitResults[fitManagerName_leg2_dR_all.Data()];
      TString outputFileName_leg2_dR_all = Form("plots/fitTauDecayKinePlots_%s_leg2_dR_all.eps", decayMode_string.Data());
      fitManager_leg2_dR_all->runPrefit(inputFileName_histograms, 
					"VisInvisDeltaRLabTimesPt", "all", outputFileName_leg2_dR_all);
      fitManager_leg2_dR_all->runFit(dataset_leg2_all, outputFileName_leg2_dR_all);
      
      TString fitManagerName_leg2_angle_all = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg2", "angle", "all");
      if ( fitResults.find(fitManagerName_leg2_angle_all.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg2_angle_all.Data()] = 
	  new fitManager(&leg2Energy, &leg2VisInvisAngleLab, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg2", "angle", "all"), momBinning);
      }
      fitManager* fitManager_leg2_angle_all = fitResults[fitManagerName_leg2_angle_all.Data()];
      TString outputFileName_leg2_angle_all = 
	Form("plots/fitTauDecayKinePlots_%s_leg2_angle_all.eps", decayMode_string.Data());
      fitManager_leg2_angle_all->runPrefit(inputFileName_histograms, 
					   "VisInvisAngleLabTimesEnergy", "all", outputFileName_leg2_angle_all);
      fitManager_leg2_angle_all->runFit(dataset_leg2_all, outputFileName_leg2_angle_all);
    }
    
    /*****************************************************************************
     ********   Run fit for leg2 with cuts on visible decay products applied  ****
     *****************************************************************************/

    if ( runSelected ) {
      RooAbsData* dataset_leg2_selected = dataset_selected->reduce(leg2DecayModeSelection.Data());
      
      TString fitManagerName_leg2_dR_selected = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg2", "dR", "selected");
      if ( fitResults.find(fitManagerName_leg2_dR_selected.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg2_dR_selected.Data()] = 
	  new fitManager(&leg2Pt, &leg2VisInvisDeltaRLab, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg2", "dR", "selected2"), momBinning);
      }
      fitManager* fitManager_leg2_dR_selected = fitResults[fitManagerName_leg2_dR_selected.Data()];
      TString outputFileName_leg2_dR_selected = Form("plots/fitTauDecayKinePlots_%s_leg2_dR_selected.eps", decayMode_string.Data());
      fitManager_leg2_dR_selected->runPrefit(inputFileName_histograms, 
					     "VisInvisDeltaRLabTimesPt", "selected2", outputFileName_leg2_dR_selected);
      fitManager_leg2_dR_selected->runFit(dataset_leg2_selected, outputFileName_leg2_dR_selected);
      
      TString fitManagerName_leg2_angle_selected = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg2", "angle", "selected");
      if ( fitResults.find(fitManagerName_leg2_angle_selected.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg2_angle_selected.Data()] = 
	  new fitManager(&leg2Energy, &leg2VisInvisAngleLab, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg2", "angle", "selected2"), momBinning);
      }
      fitManager* fitManager_leg2_angle_selected = fitResults[fitManagerName_leg2_angle_selected.Data()];
      TString outputFileName_leg2_angle_selected = 
	Form("plots/fitTauDecayKinePlots_%s_leg2_angle_selected.eps", decayMode_string.Data());
      fitManager_leg2_angle_selected->runPrefit(inputFileName_histograms, 
						"VisInvisAngleLabTimesEnergy", "selected2", outputFileName_leg2_angle_selected);
      fitManager_leg2_angle_selected->runFit(dataset_leg2_selected, outputFileName_leg2_angle_selected);
    }

//--- write results of prefit to ROOT file
    TString outputFileName_prefit = Form("fitTauDecayKinePlots_%s%s.root", decayMode_string.Data(), selection_string.Data());
    TFile* outputFile_prefit = new TFile(outputFileName_prefit.Data(), "RECREATE");
    TDirectory* dirVisInvisDeltaRLab = outputFile_prefit->mkdir("VisInvisDeltaRLabTimesPt");  
    TDirectory* dirVisInvisDeltaRLab_all = dirVisInvisDeltaRLab->mkdir("all");
    TDirectory* dirVisInvisDeltaRLab_selected = dirVisInvisDeltaRLab->mkdir("selected");
    TDirectory* dirVisInvisAngleLab = outputFile_prefit->mkdir("VisInvisAngleLabTimesEnergy");
    TDirectory* dirVisInvisAngleLab_all = dirVisInvisAngleLab->mkdir("all");
    TDirectory* dirVisInvisAngleLab_selected = dirVisInvisAngleLab->mkdir("selected");
    for ( std::map<std::string, fitManager*>::iterator prefitResult = fitResults.begin();
	  prefitResult != fitResults.end(); ++prefitResult ) {
      
      if ( !prefitResult->second ) continue;
      
      if      ( prefitResult->second->label_.Contains("_dR_")       && 
		prefitResult->second->label_.Contains("_all_")      ) dirVisInvisDeltaRLab_all->cd();
      else if ( prefitResult->second->label_.Contains("_angle_")    && 
		prefitResult->second->label_.Contains("_all_")      ) dirVisInvisAngleLab_all->cd();
      else if ( prefitResult->second->label_.Contains("_dR_")       && 
		prefitResult->second->label_.Contains("_selected_") ) dirVisInvisDeltaRLab_selected->cd();
      else if ( prefitResult->second->label_.Contains("_angle_")    && 
		prefitResult->second->label_.Contains("_selected_") ) dirVisInvisAngleLab_selected->cd();
      else assert(0);
      
      prefitResult->second->prefitResultLandauMP_->Write();
      prefitResult->second->prefitResultLandauWidth_->Write();
      prefitResult->second->prefitResultLandauScaleFactor_->Write();
      prefitResult->second->prefitResultGaussianMean_->Write();
      prefitResult->second->prefitResultGaussianSigma_->Write();
      prefitResult->second->prefitResultMix_->Write();
    }
    delete outputFile_prefit;
    
    TString outputFileName_fit = Form("fitTauDecayKinePlots_%s%s.txt", decayMode_string.Data(), selection_string.Data());
    ostream* outputFile_fit = new std::ofstream(outputFileName_fit.Data(), std::ios::out);
    for ( std::map<std::string, fitManager*>::iterator fitResult = fitResults.begin();
	  fitResult != fitResults.end(); ++fitResult ) {
      fitResult->second->writeFitResults(*outputFile_fit, fitResult->second->label_);
    }
    delete outputFile_fit;
  }

  delete dataTree;
}
