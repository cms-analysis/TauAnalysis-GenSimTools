
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

#include <Math/PdfFuncMathCore.h>

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

void pdfTimingTest(RooAbsPdf* pdf, RooRealVar* mom, RooRealVar* sepTimesMom)
{
  TRandom3 rnd;

  TBenchmark benchmark;
  benchmark.Start("pdfTimingTest");

  double pdfValueSum = 0.;

  const unsigned numCalls = 100000;
  for ( unsigned iCall = 0; iCall < numCalls; ++iCall ) {
    double momValue         = rnd.Uniform(20., 200.);
    double sepTimesMomValue = rnd.Uniform( 0.,  12.);

    mom->setVal(momValue);
    sepTimesMom->setVal(sepTimesMomValue);

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

void enlargeHistogramErrors(TH1* histogram)
{
  unsigned numBins = histogram->GetNbinsX();
  for ( unsigned iBin = 1; iBin <= numBins; ++iBin ) {
    double binContent = histogram->GetBinContent(iBin);
    double binError = histogram->GetBinError(iBin);

    double enlargedBinError = 0.1*binContent;

    histogram->SetBinError(iBin, TMath::Sqrt(binError*binError + enlargedBinError*enlargedBinError));
  }
}

double compHistogramRMSinRange(TH1* histogram, double xRef, int firstBin, int lastBin)
{
  double rms = 0.;

  double sumBinContents = 0.;

  for ( int iBin = firstBin; iBin <= lastBin; ++iBin ) {
    double binContent = histogram->GetBinContent(iBin);
    double x = histogram->GetBinCenter(iBin);

    double dx = x - xRef;
    rms += binContent*dx*dx;
     
    sumBinContents += binContent;
  }

  if ( sumBinContents != 0. ) rms = TMath::Sqrt(rms/sumBinContents);

  return rms;
}

double compHistogramRMSltMax(TH1* histogram)
{
  int binMax = histogram->GetMaximumBin();
  
  double xMax = histogram->GetBinCenter(binMax);

  return compHistogramRMSinRange(histogram, xMax, 1, binMax);
}

double compHistogramRMSgtMax(TH1* histogram)
{
  int binMax = histogram->GetMaximumBin();
  
  double xMax = histogram->GetBinCenter(binMax);

  int numBins = histogram->GetNbinsX();

  return compHistogramRMSinRange(histogram, xMax, binMax, numBins);
}

double compFallingEdgePos(TH1* histogram, double& errEdgePosRight, double& errEdgePosLeft)
{
  std::cout << "<compFallingEdgePos>:" << std::endl;

  int numBins = histogram->GetNbinsX();

  std::vector<double> window5derrivatives(numBins);

  for ( int iBin = 3; iBin <= (numBins - 2); ++iBin ) {
    double diffBinContent = histogram->GetBinContent(iBin + 2) - histogram->GetBinContent(iBin - 2);
    double diffBinCenter  = histogram->GetBinCenter(iBin + 2) - histogram->GetBinCenter(iBin - 2);

    window5derrivatives[iBin] = diffBinContent/diffBinCenter;
  }

  double minWindow5derrivative = 0;
  int binFallingEdge = -1.;
  double fallindEdgePos = -1.;
  for ( int iBin = 3; iBin <= (numBins - 2); ++iBin ) {
    if ( window5derrivatives[iBin] < minWindow5derrivative ) {
      binFallingEdge = iBin;
      fallindEdgePos = histogram->GetBinCenter(iBin);
      minWindow5derrivative = window5derrivatives[iBin];
    }
  }

  binFallingEdge -= 2; // CV: "phenomenological" correction...

  if ( binFallingEdge > 0 ) {
    errEdgePosRight = TMath::Abs(histogram->GetBinContent(binFallingEdge)/minWindow5derrivative);
    std::cout << " errEdgePosRight = " << errEdgePosRight << std::endl;
    
    double histogramMax_y = histogram->GetBinContent(histogram->GetMaximumBin());      
    errEdgePosLeft  = TMath::Abs((histogramMax_y - histogram->GetBinContent(binFallingEdge))/minWindow5derrivative);
    std::cout << " errEdgePosLeft = " << errEdgePosLeft << std::endl;
  }

  return fallindEdgePos;
}

double compFCN(TH1* histogram, int iBin0, int iBin1, double mom, double mp, double width)
{
  //std::cout << "<compFCN>:" << std::endl;

  double retVal = 0.;
  
  double histogram_sum = 0.;
  double landau_integral = 0.;
  for ( int iBin = iBin0; iBin <= iBin1; ++iBin ) {
    histogram_sum += histogram->GetBinContent(iBin);
    landau_integral += ::ROOT::Math::landau_pdf(mom*(histogram->GetBinCenter(iBin) - mp)/width);
  }

  //std::cout << " histogram_sum = " << histogram_sum << std::endl;
  //std::cout << " landau_integral = " << landau_integral << std::endl;
  
  double landau_norm = (histogram_sum/landau_integral);
  //std::cout << " landau_norm = " << landau_norm << std::endl;
  
  for ( int iBin = iBin0; iBin <= iBin1; ++iBin ) {
    double binCenter = histogram->GetBinCenter(iBin);
    
    double binContent = histogram->GetBinContent(iBin);
    //std::cout << "binContent = " << binContent << std::endl; 
    double binError = histogram->GetBinError(iBin);
    
    double landau = landau_norm*::ROOT::Math::landau_pdf(mom*(binCenter - mp)/width);
    //std::cout << "landau = " << landau << std::endl; 
    
    double pull = binContent - landau;
    if ( binError > 0. ) pull /= binError;
    //std::cout << "binCenter = " << binCenter << ": pull = " << pull << std::endl;
    
    retVal += pull*pull;
  }

  //std::cout << "--> returning retVal = " << retVal << std::endl;

  return retVal;
}

void setRealVar_Value_Range(RooRealVar* x, double value, double xMin, double xMax)
{
  bool isConstant = x->isConstant();
  x->setConstant(false);
  x->setRange(-1.e+6, +1.e+6);
  x->setVal(value);
  x->setRange(xMin, xMax);
  x->setConstant(isConstant);
}

//
//-------------------------------------------------------------------------------
//

TString convertTFormulaToRooFitSyntax(const TString& formula)
{
  std::cout << "<convertTFormulaToRooFitSyntax>:" << std::endl;
  std::cout << " formula(TFormula syntax) = " << formula.Data() << std::endl;

  TString retVal = formula;

  unsigned numParameter = 0;
  for ( unsigned iParameter = 0; iParameter < 10; iParameter++ ) {
    TString parameterNameTFormula = Form("[%u]", iParameter);
    TString parameterNameRooFit = Form("@%u", iParameter);

    if ( retVal.Contains(parameterNameTFormula) ) {
      retVal = retVal.ReplaceAll(parameterNameTFormula, parameterNameRooFit);
      numParameter = iParameter + 1;
    }
  }

  TString xNameTFormula = "x";
  TString xNameRooFit = Form("@%u", numParameter);

  int index = retVal.Index(xNameTFormula.Data());
  while ( index != kNPOS ) {
//--- skip replacing 'x' in case it is part of a character sequence
//   (representing a function call)
    if ( (index == 0                     || !TString(retVal[index - 1]).IsAlpha()) &&
	 (index == (retVal.Length() - 1) || !TString(retVal[index + 1]).IsAlpha()) )
      retVal = retVal.Replace(index, xNameTFormula.Length(), xNameRooFit.Data());
    index = retVal.Index(xNameTFormula.Data(), index + xNameRooFit.Length());
  }

  std::cout << "--> formula(RooFit syntax) = " << retVal.Data() << std::endl;

  return retVal;
}

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
  fitManager(RooRealVar* mom, RooRealVar* sep, RooRealVar* sepTimesMom, 
	     int decayMode, const TString& label, const TArrayD& momBinning)
    : mom_(mom),
      sep_(sep),
      sepTimesMom_(sepTimesMom),
      decayMode_(decayMode),
      label_(label),
      momBinning_(momBinning),
      prefitParamGMean_(momBinning.GetSize() - 1),
      prefitParamGSigma_(momBinning.GetSize() - 1),
      prefitParamAlpha_(momBinning.GetSize() - 1),
      prefitParamSlope_(momBinning.GetSize() - 1),
      prefitParamOffset_(momBinning.GetSize() - 1),
      prefitParamC_(momBinning.GetSize() - 1),
      prefitParamMP1_(momBinning.GetSize() - 1),
      prefitParamWidth1_(momBinning.GetSize() - 1),
      prefitParamMP2_(momBinning.GetSize() - 1),
      prefitParamWidth2_(momBinning.GetSize() - 1),
      prefitParamX0_(momBinning.GetSize() - 1),
      prefitParamDeltaX1_(momBinning.GetSize() - 1),
      prefitModel_(momBinning.GetSize() - 1),
      fitParamCoeffGMean_(5),
      fitParamGMean_(0),
      fitParamCoeffGSigma_(5),      
      fitParamGSigma_(0),
      fitParamCoeffAlpha_(2),      
      fitParamAlpha_(0),
      fitParamCoeffSlope_(2),
      fitParamSlope_(0),
      fitParamCoeffOffset_(2),
      fitParamOffset_(0),
      fitParamCoeffC_(4),
      fitParamC_(0),
      fitParamCoeffMP1_(5),
      fitParamMP1_(0),
      fitParamCoeffWidth1_(2),
      fitParamWidth1_(0),
      fitParamCoeffMP2_(3),
      fitParamMP2_(0),
      fitParamCoeffWidth2_(2),
      fitParamWidth2_(0),
      fitParamCoeffX0_(2),
      fitParamX0_(0),
      fitParamCoeffDeltaX1_(3),
      fitParamDeltaX1_(0),
      fitModel_(0)
  {
    std::cout << "<fitManager::fitManager>:" << std::endl;
    std::cout << " mom = " << mom_->GetName() << std::endl;
    std::cout << " dR = " << sep_->GetName() << std::endl;
    std::cout << " decayMode = " << getDecayMode_string(decayMode_) << std::endl;

    prefitResultGMean_ = new TGraphErrors();
    prefitResultGSigma_ = new TGraphErrors();
    prefitResultAlpha_ = new TGraphErrors();
    prefitResultSlope_ = new TGraphErrors();
    prefitResultOffset_ = new TGraphErrors();
    prefitResultC_ = new TGraphErrors();
    prefitResultMP1_ = new TGraphErrors();
    prefitResultWidth1_ = new TGraphErrors();
    prefitResultMP2_ = new TGraphErrors();
    prefitResultWidth2_ = new TGraphErrors();
    prefitResultX0_ = new TGraphErrors();
    prefitResultDeltaX1_ = new TGraphErrors();
  
//--- parametrize pt/energy dependence of GMean fitParameter by (orthogonal) Chebyshev polynomials
//   ( cf. http://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html )
    momParametrizationFormulaGMean_   = 
      "TMath::Abs((1./(x*x*x*x))*([0] + [1]*x + [2]*(2.0*x*x - 1.0) + [3]*(4.0*x*x*x - 3.0*x) + [4]*(8.0*x*x*x*x - 8.0*x*x + 1.0)))";
    momParametrizationFormulaGSigma_  = "([0] + [1]*x)*(1.0 + [2]*TMath::TanH([3] + [4]*x))";
    momParametrizationFormulaAlpha_   = "[0] + [1]*x";
    momParametrizationFormulaSlope_   = "[0] + [1]*x";
    momParametrizationFormulaOffset_  = "[0] + [1]*x";
    momParametrizationFormulaC_       = 
      "0.5*(1.0 + TMath::Cos([0])) + 0.125*(1.0 - TMath::Cos([0]))*(1.0 + TMath::Cos([1]))*(1.0 + TMath::TanH([2] + [3]*x))";
//--- parametrize pt/energy dependence of MP1 fitParameter by (orthogonal) Chebyshev polynomials
//   ( cf. http://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html )
    momParametrizationFormulaMP1_     =
      "TMath::Abs((1./(x*x*x*x))*([0] + [1]*x + [2]*(2.0*x*x - 1.0) + [3]*(4.0*x*x*x - 3.0*x) + [4]*(8.0*x*x*x*x - 8.0*x*x + 1.0)))";
    momParametrizationFormulaWidth1_  = "TMath::Min(1.e-4, [0] + [1]*x)";
    momParametrizationFormulaMP2_     = "[0]*TMath::Min(1., [1]*(x - [2]))";
    momParametrizationFormulaWidth2_  = "TMath::Min(5.e-2, [0] + [1]*x)";
    momParametrizationFormulaX0_      = "[0] + [1]*x";
    momParametrizationFormulaDeltaX1_ = "[0]*TMath::Min(1., [1]*(x - [2]))";
    
    momParamGMean_ = bookTF1("GMean", momParametrizationFormulaGMean_);
    momParamGMean_->SetParameter(   0,      2.9        );
    momParamGMean_->SetParameter(   1,      0.         );
    momParamGMean_->SetParameter(   2,      0.2        );
    momParamGMean_->SetParameter(   3, -60.*0.05       );
    momParamGMean_->SetParameter(   4,      0.05       );
    momParamGSigma_ = bookTF1("GSigma", momParametrizationFormulaGSigma_);
    momParamGSigma_->SetParameter(  0,      1.         );
    momParamGSigma_->SetParameter(  1,      0.         );
    momParamGSigma_->SetParameter(  2,      0.2        );
    momParamGSigma_->SetParameter(  3, -60.*0.05       );
    momParamGSigma_->SetParameter(  4,      0.05       );
    momParamAlpha_ = bookTF1("Alpha", momParametrizationFormulaAlpha_);
    momParamSlope_ = bookTF1("Slope", momParametrizationFormulaSlope_);
    momParamOffset_ = bookTF1("Offset", momParametrizationFormulaOffset_);
    momParamC_ = bookTF1("C", momParametrizationFormulaC_);
    momParamC_->SetParameter(       0, 0.5*TMath::Pi() );
    momParamC_->SetParameter(       1, 0.5*TMath::Pi() );
    momParamC_->SetParameter(       2, -60.*0.05       );
    momParamC_->SetParameter(       3,      0.05       );
    momParamMP1_ = bookTF1("MP1", momParametrizationFormulaMP1_);
    momParamWidth1_ = bookTF1("Width1", momParametrizationFormulaWidth1_);
    momParamMP2_ = bookTF1("MP2", momParametrizationFormulaMP2_);
    momParamMP2_->SetParameter(     0,      9.5        );
    momParamMP2_->SetParameter(     1,      0.01       );
    momParamMP2_->SetParameter(     2,     10.         );
    momParamWidth2_ = bookTF1("Width2", momParametrizationFormulaWidth2_);
    momParamX0_ = bookTF1("X0", momParametrizationFormulaX0_);
    momParamDeltaX1_ = bookTF1("DeltaX1", momParametrizationFormulaDeltaX1_);
    momParamDeltaX1_->SetParameter( 0,      8.5        );
    momParamDeltaX1_->SetParameter( 1,      0.01       );
    momParamDeltaX1_->SetParameter( 2,     20.         );
  }

  ~fitManager()
  {
    clearCollection(prefitParamGMean_);
    clearCollection(prefitParamGSigma_);
    clearCollection(prefitParamAlpha_);
    clearCollection(prefitParamSlope_);
    clearCollection(prefitParamOffset_);
    clearCollection(prefitParamC_);
    clearCollection(prefitParamMP1_);
    clearCollection(prefitParamWidth1_);
    clearCollection(prefitParamMP2_);
    clearCollection(prefitParamWidth2_);
    clearCollection(prefitParamX0_);
    clearCollection(prefitParamDeltaX1_);
    clearCollection(prefitModel_);

    delete prefitResultGMean_;
    delete prefitResultGSigma_;
    delete prefitResultAlpha_;
    delete prefitResultSlope_;
    delete prefitResultOffset_;
    delete prefitResultC_;
    delete prefitResultMP1_;
    delete prefitResultWidth1_;
    delete prefitResultMP2_;
    delete prefitResultWidth2_;
    delete prefitResultX0_;
    delete prefitResultDeltaX1_;

    delete momParamGMean_;
    delete momParamGSigma_;
    delete momParamAlpha_;
    delete momParamSlope_;
    delete momParamOffset_;
    delete momParamC_;
    delete momParamMP1_;
    delete momParamWidth1_;
    delete momParamMP2_;
    delete momParamWidth2_;
    delete momParamX0_;
    delete momParamDeltaX1_;
    
    clearCollection(fitParamCoeffGMean_);
    delete fitParamGMean_;
    clearCollection(fitParamCoeffGSigma_);
    delete fitParamGSigma_;
    clearCollection(fitParamCoeffAlpha_);
    delete fitParamAlpha_;
    clearCollection(fitParamCoeffSlope_);
    delete fitParamSlope_;
    clearCollection(fitParamCoeffOffset_);
    delete fitParamOffset_;
    clearCollection(fitParamCoeffC_);
    delete fitParamC_;
    clearCollection(fitParamCoeffMP1_);
    delete fitParamMP1_;
    clearCollection(fitParamCoeffWidth1_);
    delete fitParamWidth1_;
    clearCollection(fitParamCoeffMP2_);
    delete fitParamMP2_;
    clearCollection(fitParamCoeffWidth2_);
    delete fitParamWidth2_;
    clearCollection(fitParamCoeffX0_);
    delete fitParamX0_;
    clearCollection(fitParamCoeffDeltaX1_);
    delete fitParamDeltaX1_;
    delete fitModel_;
  }

  TF1* bookTF1(const TString& fitParameter, const TString& formula)
  {
    std::cout << "<fitManager::bookTF1>:" << std::endl;

    TString decayMode_string = getDecayMode_string(decayMode_);
    TString tf1Name = Form("%s_%s_%s_%s", fitParameter.Data(), sep_->GetName(), decayMode_string.Data(), label_.Data());
    double momMin = momBinning_[0];
    //std::cout << "momMin = " << momMin << std::endl;
    double momMax = momBinning_[momBinning_.GetSize() - 1];
    //std::cout << "momMax = " << momMax << std::endl;
    TF1* tf1 = new TF1(tf1Name.Data(), formula, momMin, momMax);
    return tf1;
  }

  void fitTF1(TF1* tf1, TGraph* graph, std::vector<double>& fitParamValues, std::vector<double>& fitParamErrors)
  {
    std::cout << "<fitManager::fitTF1>:" << std::endl;

    int numFitParameter = tf1->GetNpar();
    std::cout << " numFitParameter = " << numFitParameter << std::endl;

    bool isFitted = false;
    if ( graph->GetN() >= (2*numFitParameter) ) {
      graph->Fit(tf1, "W0");
      isFitted = true;
    } else {
      std::cout << "Number of points in Graph  = " << graph->GetN() << ","
		<< " insufficient for fitting Formula = " << tf1->GetTitle() 
		<< " with " << tf1->GetNumberFreeParameters() << " free Parameters --> skipping TF1 fit !!" << std::endl;
    }
    
    // CV: if TGraph has been fitted, 
    //     set starting-value for final fit to fitted TF1 parameters,
    //     else set starting-value to initial TF1 parameters
    fitParamValues.resize(numFitParameter);
    fitParamErrors.resize(numFitParameter);

    for ( int iParameter = 0; iParameter < numFitParameter; ++iParameter ) {
      fitParamValues[iParameter] = tf1->GetParameter(iParameter);
      fitParamErrors[iParameter] = tf1->GetParError(iParameter);
    }

    if ( isFitted ) {
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
      
      TString outputFileName = TString("plots/fitTF1").Append(tf1->GetName()).Append(".eps");
      
      canvas->SaveAs(outputFileName.Data());
      canvas->SaveAs(outputFileName.ReplaceAll(".eps", ".root").Data());
      
      delete canvas;
    }
  }

  RooRealVar* buildPrefitParameter(const TString& name, const TString& momBinName, double startValue, double min, double max)
  {
    std::cout << "<fitManager::buildPrefitParameter>:" << std::endl;

    TString decayMode_string = getDecayMode_string(decayMode_);
    TString fitParameterName = 
      Form("%s_%s_%s_%s_%s", name.Data(), momBinName.Data(), sep_->GetName(), decayMode_string.Data(), label_.Data());
    RooRealVar* fitParameter = new RooRealVar(fitParameterName.Data(), fitParameterName.Data(), startValue, min, max);
    return fitParameter;
  }

  RooConstVar* buildConstPrefitParameter(const TString& name, const TString& momBinName, double value)
  {
    std::cout << "<buildConstPrefitParameter>:" << std::endl;

    TString decayMode_string = getDecayMode_string(decayMode_);
    TString fitParameterName = 
      Form("%s_%s_%s_%s_%s", name.Data(), momBinName.Data(), sep_->GetName(), decayMode_string.Data(), label_.Data());
    RooConstVar* fitParameter = new RooConstVar(fitParameterName.Data(), fitParameterName.Data(), value);
    return fitParameter;
  }

  void printPrefitParameter(const TString& fitParameterName, RooRealVar* fitParameter)
  {
    std::cout << fitParameterName.Data() << ": " << fitParameter->getVal() << " +/- " << fitParameter->getError() << std::endl;
  }

  void printAllPrefitParameter(int momBin)
  {
    printPrefitParameter("gmean", prefitParamGMean_[momBin]);
    printPrefitParameter("gsigma", prefitParamGSigma_[momBin]);
    printPrefitParameter("alpha", prefitParamAlpha_[momBin]);
    printPrefitParameter("slope", prefitParamSlope_[momBin]);
    printPrefitParameter("offset", prefitParamOffset_[momBin]);
    printPrefitParameter("C", prefitParamC_[momBin]);
    printPrefitParameter("mp1", prefitParamMP1_[momBin]);
    printPrefitParameter("width1", prefitParamWidth1_[momBin]);
    printPrefitParameter("mp2", prefitParamMP2_[momBin]);
    printPrefitParameter("width2", prefitParamWidth2_[momBin]);
    printPrefitParameter("x0", prefitParamX0_[momBin]);
    printPrefitParameter("dx1", prefitParamDeltaX1_[momBin]);
  }

  void runPrefit(const TString& inputFileName, const TString& inputDirName, const TString& label, const TString& outputFileName)
  {
    std::cout << "<fitManager::runPrefit>:" << std::endl;
    std::cout << " inputFileName = " << inputFileName.Data() << std::endl;
    std::cout << " inputDirName = " << inputDirName.Data() << std::endl;
    
    printTimeStamp("<fitManager::runPrefit>");

    TString decayMode_string = getDecayMode_string(decayMode_);

    RooRealVar* sepTimesMom_prefit = new RooRealVar("sepTimesMom_prefit", "sepTimesMom_prefit", 0., 12.);
    sepTimesMom_prefit->setBins(120);

    TString inputDirName_full = Form("%s/%s/%s", inputDirName.Data(), label.Data(), decayMode_string.Data());
    std::cout << " inputDirName_full = " << inputDirName_full.Data() << std::endl;

    TFile* inputFile = new TFile(inputFileName.Data());

    TCanvas* canvas = new TCanvas("canvas", "canvas", 1, 1, 800, 600);
    canvas->SetFillColor(10);
    canvas->SetBorderSize(2);
  
    for ( int iMomBin = 0; iMomBin < (momBinning_.GetSize() - 1); ++iMomBin ) {
      double momMin = momBinning_[iMomBin];
      double momMax = momBinning_[iMomBin + 1];
    
      TString momName_string = "Mom";
      if      ( inputDirName.Contains("DeltaR") ) momName_string = "Pt";
      else if ( inputDirName.Contains("Angle")  ) momName_string = "Energy";

      TString histogramName = 
	Form("%s_%s%s%2.0fto%2.0f", decayMode_string.Data(), label.Data(), momName_string.Data(), momMin, momMax);
      TH1* histogram = (TH1*)inputFile->Get(TString(inputDirName_full.Data()).Append("/").Append(histogramName.Data()));      
      std::cout << " histogramName = " << histogramName.Data() << ": histogram = " << histogram << std::endl;
      enlargeHistogramErrors(histogram);

      if ( histogram->GetEntries() < 500 ) continue;

      TString datahistName = TString(histogramName).Append("_datahist");
      RooDataHist* datahist = new RooDataHist(datahistName.Data(), datahistName.Data(), RooArgList(*sepTimesMom_prefit), histogram);
      std::cout << " datahist = " << datahist << std::endl;

      TString momBinName = Form("%s%2.0fto%2.0f", momName_string.Data(), momMin, momMax);

      double histogramMean = histogram->GetMean();
      //if ( histogramMean > 0.15 ) histogramMean = 0.15;
      std::cout << " histogramMean = " << histogramMean << std::endl;
      double histogramMax_x = histogram->GetBinCenter(histogram->GetMaximumBin());
      //double histogramMax_y = histogram->GetBinContent(histogram->GetMaximumBin());
      //double risingEdge_slope = histogramMax_y/histogramMax_x;
      //if ( histogramMax_x > 0.15 ) histogramMax = 0.15;
      std::cout << " histogramMax_x = " << histogramMax_x << std::endl;
      double errEdgePosRight = 0.;
      double errEdgePosLeft  = 0.;
      double fallingEdge_position = compFallingEdgePos(histogram, errEdgePosRight, errEdgePosLeft);
      std::cout << " fallingEdge_position = " << fallingEdge_position << std::endl;    
      double histogramRMS = histogram->GetRMS();
      //if ( histogramRMS  > 0.15 ) histogramRMS  = 0.15;
      std::cout << " histogramRMS = " << histogramRMS << std::endl;
      double histogramRMSltMax = compHistogramRMSltMax(histogram);
      std::cout << " histogramRMSltMax = " << histogramRMSltMax << std::endl;
      double histogramRMSgtMax = compHistogramRMSgtMax(histogram);
      std::cout << " histogramRMSgtMax = " << histogramRMSgtMax << std::endl;

      const double epsilon = 1.e-3;

      prefitParamGMean_[iMomBin] = buildPrefitParameter("gmean", momBinName, histogramMax_x, 0.5*histogramMax_x, 1.5*histogramMax_x);
      double gsigmaMin = 0.25*TMath::Min(histogramRMSltMax, histogramRMS);
      double gsigmaMax = 2.0*TMath::Max(histogramRMSltMax, histogramRMS);
      prefitParamGSigma_[iMomBin] = buildPrefitParameter("gsigma", momBinName, histogramRMS, gsigmaMin, gsigmaMax);
      //prefitParamAlpha_[iMomBin] = buildPrefitParameter("alpha", momBinName, 0., -1.e+2, +1.e-0);
      prefitParamAlpha_[iMomBin] = buildPrefitParameter("alpha", momBinName, 0., -1.e-9, +1.e-9);
      prefitParamAlpha_[iMomBin]->setConstant();
      double slopeMin = 1.e-2/histogramMax_x;
      double slopeMax = 1.e+2/histogramMax_x;
      prefitParamSlope_[iMomBin] = buildPrefitParameter("slope", momBinName, 1./histogramMax_x, slopeMin, slopeMax);
      //prefitParamSlope_[iMomBin] = buildPrefitParameter("slope", momBinName, 0., -1.e-9, +1.e-9);
      //prefitParamSlope_[iMomBin]->setConstant();
      //prefitParamOffset_[iMomBin] = buildPrefitParameter("offset", momBinName, 1.e-9, 0., +0.1);
      bool offset_isConstant = true;
      if ( decayMode_string == "OneProng0Pi0" ) {
	prefitParamOffset_[iMomBin] = buildPrefitParameter("offset", momBinName, 0., -1.e+2, +1.e+2);
	offset_isConstant = false;
      } else {
	prefitParamOffset_[iMomBin] = buildPrefitParameter("offset", momBinName, 0., -1.e-9, +1.e-9);
	prefitParamOffset_[iMomBin]->setConstant();
      }
      prefitParamC_[iMomBin] = buildPrefitParameter("C", momBinName, 0.25, 0., 1.);  
      //prefitParamC_[iMomBin] = buildPrefitParameter("C", momBinName, 1., 1. - epsilon, 1. + epsilon);
      //prefitParamC_[iMomBin]->setConstant();
      prefitParamMP1_[iMomBin] = buildPrefitParameter("mp1", momBinName, histogramMax_x, 0.8*histogramMax_x, fallingEdge_position);
      prefitParamWidth1_[iMomBin] = buildPrefitParameter("width1", momBinName, 0.04*histogramRMSgtMax, 1.e-4, 2.0*histogramRMSgtMax);
      //prefitParamWidth1_[iMomBin] = buildPrefitParameter("width1", momBinName, 1.e-3, 1.e-3 - epsilon, 1.e-3 + epsilon); 
      //prefitParamWidth1_[iMomBin]->setConstant();
      prefitParamMP2_[iMomBin] = buildPrefitParameter("mp2", momBinName, histogramMax_x + 2., fallingEdge_position + 2., 12.);
      prefitParamWidth2_[iMomBin] = buildPrefitParameter("width2", momBinName, 5.e-2, 5.e-3, 2.e-1);
      double x0Min = 0.9*fallingEdge_position - errEdgePosLeft;
      double x0Max = 1.1*fallingEdge_position + errEdgePosRight;
      prefitParamX0_[iMomBin] = buildPrefitParameter("x0", momBinName, fallingEdge_position, x0Min, x0Max);
      prefitParamDeltaX1_[iMomBin] = buildPrefitParameter("dx1", momBinName, 12. - epsilon, 0., 12.);

      TString modelName = Form("pdf_%s_%s_%s_%s", decayMode_string.Data(), momBinName.Data(), sep_->GetName(), label_.Data());
      prefitModel_[iMomBin] = 
	new TauDecayKinePdf(modelName.Data(), modelName.Data(), 
			    *sepTimesMom_prefit, 
			    *prefitParamGMean_[iMomBin], *prefitParamGSigma_[iMomBin], *prefitParamAlpha_[iMomBin], 
			    *prefitParamSlope_[iMomBin], *prefitParamOffset_[iMomBin], *prefitParamC_[iMomBin],
			    *prefitParamMP1_[iMomBin], *prefitParamWidth1_[iMomBin], 
			    *prefitParamMP2_[iMomBin], *prefitParamWidth2_[iMomBin], 
			    *prefitParamX0_[iMomBin], *prefitParamDeltaX1_[iMomBin]);

      prefitModel_[iMomBin]->disableAnalyticIntegration();

      RooConstVar* gmeanConstraint_value = 
        new RooConstVar("gmeanConstraint_value", "gmeanConstraint_value", histogramMax_x);
      RooConstVar* gmeanConstraint_sigma =
        new RooConstVar("gmeanConstraint_sigma", "gmeanConstraint_sigma", 0.25*histogramMax_x);
      RooGaussian* gmeanConstraint_pdf =
        new RooGaussian("gmeanConstraint_pdf", "gmeanConstraint_pdf",
			*prefitParamGMean_[iMomBin], *gmeanConstraint_value, *gmeanConstraint_sigma);

      RooConstVar* slopeConstraint_value = 
        new RooConstVar("slopeConstraint_value", "slopeConstraint_value", 1./histogramMax_x);
      RooConstVar* slopeConstraint_sigma =
        new RooConstVar("slopeConstraint_sigma", "slopeConstraint_sigma", 10./histogramMax_x);
      RooGaussian* slopeConstraint_pdf =
        new RooGaussian("slopeConstraint_pdf", "slopeConstraint_pdf",
			*prefitParamSlope_[iMomBin], *slopeConstraint_value, *slopeConstraint_sigma);

      RooConstVar* x0Constraint_value =
        new RooConstVar("x0Constraint_value", "x0Constraint_value", histogramMax_x);
      RooConstVar* x0Constraint_sigma =
        new RooConstVar("x0Constraint_sigma", "x0Constraint_sigma", 0.5*histogramMax_x);
      RooGaussian* x0Constraint_pdf =
        new RooGaussian("x0Constraint_pdf", "x0Constraint_pdf",
			*prefitParamX0_[iMomBin], *x0Constraint_value, *x0Constraint_sigma);

      RooConstVar* dx1Constraint_value =
        new RooConstVar("dx1Constraint_value", "dx1Constraint_value", 12.);
      RooConstVar* dx1Constraint_sigma =
        new RooConstVar("dx1Constraint_sigma", "dx1Constraint_sigma", 12.);
      RooGaussian* dx1Constraint_pdf =
        new RooGaussian("dx1Constraint_pdf", "dx1Constraint_pdf",
			*prefitParamDeltaX1_[iMomBin], *dx1Constraint_value, *dx1Constraint_sigma);
      
      RooArgSet prefitParamConstraints(*gmeanConstraint_pdf, *slopeConstraint_pdf, *x0Constraint_pdf, *dx1Constraint_pdf);

      RooLinkedList options;
      options.Add(new RooCmdArg(RooFit::Save(true)));
      //options.Add(new RooCmdArg(RooFit::SumW2Error(true)));
      options.Add(new RooCmdArg(RooFit::PrintLevel(-1)));
      options.Add(new RooCmdArg(RooFit::PrintEvalErrors(-1)));
      options.Add(new RooCmdArg(RooFit::Warnings(-1)));
      options.Add(new RooCmdArg(RooFit::ExternalConstraints(prefitParamConstraints)));

//--- perform stand-alone fit of Gaussian/ExpGamma distribution
      std::cout << "--> fitting Gaussian/ExpGamma distribution..." << std::endl;
      RooLinkedList options_gaussian(options);
      sepTimesMom_prefit->setRange("gaussian", 0., fallingEdge_position);
      options_gaussian.Add(new RooCmdArg(RooFit::Range("gaussian")));
      prefitParamGMean_[iMomBin]->setConstant(false);
      prefitParamGSigma_[iMomBin]->setConstant(false);
      //prefitParamAlpha_[iMomBin]->setConstant(false);
      prefitParamSlope_[iMomBin]->setConstant(false);
      if ( !offset_isConstant ) prefitParamOffset_[iMomBin]->setConstant(false);
      prefitParamC_[iMomBin]->setConstant(false);
      prefitParamMP1_[iMomBin]->setConstant(true);
      prefitParamWidth1_[iMomBin]->setConstant(true);
      prefitParamMP2_[iMomBin]->setConstant(true);
      prefitParamWidth2_[iMomBin]->setConstant(true);
      prefitParamX0_[iMomBin]->setConstant(true);
      prefitParamDeltaX1_[iMomBin]->setConstant(true);
      RooFitResult* prefitResult_gaussian = prefitModel_[iMomBin]->fitTo(*datahist, options_gaussian);
      delete prefitResult_gaussian;

      printAllPrefitParameter(iMomBin);

//--- perform stand-alone fit of Landau distribution
      std::cout << "--> fitting Landau1 distribution..." << std::endl;
      RooLinkedList options_landau1(options);
      int iBin0 = histogram->FindBin(0.5*(histogramMax_x + fallingEdge_position));
      const int numBinsRef = 5;
      int iBinRef = histogram->FindBin(fallingEdge_position) + numBinsRef;
      std::cout << "iBinRef = " << iBinRef << ": position = " << histogram->GetBinCenter(iBinRef) << std::endl;
      sepTimesMom_prefit->setRange("landau1", histogram->GetBinCenter(iBin0), histogram->GetBinCenter(iBinRef));
      options_landau1.Add(new RooCmdArg(RooFit::Range("landau1")));
      prefitParamGMean_[iMomBin]->setConstant(true);
      prefitParamGSigma_[iMomBin]->setConstant(true);
      //prefitParamAlpha_[iMomBin]->setConstant(true);
      prefitParamSlope_[iMomBin]->setConstant(true);
      if ( !offset_isConstant ) prefitParamOffset_[iMomBin]->setConstant(true);
      prefitParamC_[iMomBin]->setConstant(true);
      prefitParamMP1_[iMomBin]->setConstant(false);
      prefitParamWidth1_[iMomBin]->setConstant(false);
      prefitParamMP2_[iMomBin]->setConstant(true);
      prefitParamWidth2_[iMomBin]->setConstant(true);
      prefitParamX0_[iMomBin]->setConstant(true);
      prefitParamDeltaX1_[iMomBin]->setConstant(true);
      RooFitResult* prefitResult_landau1 = prefitModel_[iMomBin]->fitTo(*datahist, options_landau1);
//--- CV: return value of RooFitResult::minNll is unreliable,
//        because it may change by large amount in case Hesse matrix has negative diagonal elements
//       --> call dedicated function in order to compute FCN value a-la Minuit
      double refFCN = compFCN(histogram, iBin0, iBinRef, 0.5*(momMax + momMin), 
			      prefitParamMP1_[iMomBin]->getVal(), prefitParamWidth1_[iMomBin]->getVal());
      //std::cout << "refFCN = " << refFCN << std::endl;
      double refFCNperDoF = refFCN/((iBinRef - iBin0) + 1 - 2);
      //std::cout << "refFCNperDoF = " << refFCNperDoF << std::endl;
      delete prefitResult_landau1;
      int iBin1 = histogram->GetNbinsX();
      //std::cout << "iBin1 = " << iBin1 << ": position = " << histogram->GetBinCenter(iBin1) << std::endl;
      int dBin = TMath::CeilNint(0.5*(histogram->GetNbinsX() - iBin0));
      //std::cout << "dBin = " << dBin << std::endl;
      double FCNperDoF = 0.;
      do {	
	sepTimesMom_prefit->setRange("landau1", histogram->GetBinCenter(iBin0), histogram->GetBinCenter(iBin1));
	RooFitResult* prefitResult_landau1 = prefitModel_[iMomBin]->fitTo(*datahist, options_landau1);
	double FCN = compFCN(histogram, iBin0, iBin1, 0.5*(momMax + momMin),
			     prefitParamMP1_[iMomBin]->getVal(), prefitParamWidth1_[iMomBin]->getVal());
	//std::cout << "FCN = " << FCN << std::endl;
	FCNperDoF = FCN/((iBin1 - iBin0) + 1 - 2);
	//std::cout << "FCNperDoF = " << FCNperDoF << std::endl;
	delete prefitResult_landau1;	
	if ( FCNperDoF > (2.*refFCNperDoF) ) iBin1 = TMath::Max(iBinRef, iBin1 - dBin);
	else iBin1 = TMath::Min(histogram->GetNbinsX(), iBin1 + dBin);
	//std::cout << "iBin1 = " << iBin1 << ": position = " << histogram->GetBinCenter(iBin1) << std::endl;
	dBin = TMath::CeilNint(0.5*dBin);
	//std::cout << "dBin = " << dBin << std::endl;
      } while ( dBin > 1 );

      if ( FCNperDoF > (2.*refFCNperDoF) ) iBin1 = TMath::Max(iBinRef, iBin1 - dBin);
      iBin1 -= 3; // CV: "phenomenological" correction...

      printAllPrefitParameter(iMomBin);
      
      double x1 = histogram->GetBinCenter(iBin1);
      if ( iBin1 > (histogram->GetNbinsX() - 5) ) x1 = 12.;
      std::cout << "x1 = " << x1 << std::endl;

      setRealVar_Value_Range(prefitParamMP2_[iMomBin], 0.8*x1, 0.6*x1, x1);
      double dx1Min = histogram->GetBinCenter(iBin1 - 2);
      double dx1Max = histogram->GetBinCenter(iBin1 + 2);
      if ( dx1Max > 12. ) dx1Max = 12.;
      setRealVar_Value_Range(prefitParamDeltaX1_[iMomBin], x1 - fallingEdge_position, dx1Min, dx1Max);

//--- perform stand-alone fit of 'Exponential' distribution
      if ( x1 < 12. ) {
	std::cout << "--> fitting Landau2 distribution..." << std::endl;
	RooLinkedList options_landau2(options);
	sepTimesMom_prefit->setRange("landau2", x1, 12.);
	options_landau2.Add(new RooCmdArg(RooFit::Range("landau2")));
	prefitParamGMean_[iMomBin]->setConstant(true);
	prefitParamGSigma_[iMomBin]->setConstant(true);
	//prefitParamAlpha_[iMomBin]->setConstant(true);	
	prefitParamSlope_[iMomBin]->setConstant(true);
	if ( !offset_isConstant ) prefitParamOffset_[iMomBin]->setConstant(true);
	prefitParamC_[iMomBin]->setConstant(true);
	prefitParamMP1_[iMomBin]->setConstant(true);
	prefitParamWidth1_[iMomBin]->setConstant(true);
	prefitParamMP2_[iMomBin]->setConstant(false);
	prefitParamWidth2_[iMomBin]->setConstant(false);
	prefitParamX0_[iMomBin]->setConstant(true);
	prefitParamDeltaX1_[iMomBin]->setConstant(true);
	RooFitResult* prefitResult_landau2 = prefitModel_[iMomBin]->fitTo(*datahist, options_landau2);
	delete prefitResult_landau2;
	
	printAllPrefitParameter(iMomBin);
      } else {
	std::cout << "--> skipping fit of Landau2 distribution..." << std::endl;
      }

//--- start combined fit of Gaussian/ExpGamma + Landau + Exponential model
      std::cout << "--> starting combined fit of Gaussian/ExpGamma + Landau1 + Landau2 model..." << std::endl;
      RooLinkedList options_combined(options);
      sepTimesMom_prefit->setRange("combined", 0., 12.);
      options_combined.Add(new RooCmdArg(RooFit::Range("combined")));
      prefitParamGMean_[iMomBin]->setConstant(false);
      prefitParamGSigma_[iMomBin]->setConstant(false);
      //prefitParamAlpha_[iMomBin]->setConstant(false);
      prefitParamSlope_[iMomBin]->setConstant(false);
      if ( !offset_isConstant ) prefitParamOffset_[iMomBin]->setConstant(false);
      prefitParamC_[iMomBin]->setConstant(false);
      prefitParamMP1_[iMomBin]->setConstant(false);
      prefitParamWidth1_[iMomBin]->setConstant(false);
      prefitParamMP2_[iMomBin]->setConstant(false);
      prefitParamWidth2_[iMomBin]->setConstant(false);
      prefitParamX0_[iMomBin]->setConstant(false);
      prefitParamDeltaX1_[iMomBin]->setConstant(false);
      RooFitResult* prefitResult = prefitModel_[iMomBin]->fitTo(*datahist, options_combined);
      std::cout << " prefit status = " << prefitResult->status() << " (converged = 0)" << std::endl;
      delete prefitResult;
  
//--- CV: map first Landau --> second Landau in case only one Landau distribution is needed to fit the TAUOLA prediction,
//        in order to make dx1 linearly increasing at low pt/energy
      if ( (prefitParamX0_[iMomBin]->getVal() + prefitParamDeltaX1_[iMomBin]->getVal()) > 11.5 &&
	   ((0.5*(momMin + momMax) < 75. && decayMode_string == "Electron_Muon")  || 
	    (0.5*(momMin + momMax) < 75. && decayMode_string == "ThreeProng0Pi0")) ) {
	std::cout << "mapping Landau1 --> Landau2." << std::endl;
	setRealVar_Value_Range(prefitParamMP2_[iMomBin], prefitParamMP1_[iMomBin]->getVal(), 
			       prefitParamMP1_[iMomBin]->getVal() - epsilon, prefitParamMP1_[iMomBin]->getVal() + epsilon);
	setRealVar_Value_Range(prefitParamWidth2_[iMomBin], prefitParamWidth1_[iMomBin]->getVal(), 
			       prefitParamWidth1_[iMomBin]->getVal() - epsilon, prefitParamWidth1_[iMomBin]->getVal() + epsilon);
	setRealVar_Value_Range(prefitParamDeltaX1_[iMomBin], 0., -epsilon, +epsilon);
      }

      printAllPrefitParameter(iMomBin);

      prefitModel_[iMomBin]->enableAnalyticIntegration();

      storePrefitResults(prefitParamGMean_[iMomBin], momMin, momMax, prefitResultGMean_);
      storePrefitResults(prefitParamGSigma_[iMomBin], momMin, momMax, prefitResultGSigma_);
      storePrefitResults(prefitParamAlpha_[iMomBin], momMin, momMax, prefitResultAlpha_);
      storePrefitResults(prefitParamSlope_[iMomBin], momMin, momMax, prefitResultSlope_);
      storePrefitResults(prefitParamOffset_[iMomBin], momMin, momMax, prefitResultOffset_);
      storePrefitResults(prefitParamC_[iMomBin], momMin, momMax, prefitResultC_);
      storePrefitResults(prefitParamMP1_[iMomBin], momMin, momMax, prefitResultMP1_);
      storePrefitResults(prefitParamWidth1_[iMomBin], momMin, momMax, prefitResultWidth1_);
      if ( (prefitParamX0_[iMomBin]->getVal() + prefitParamDeltaX1_[iMomBin]->getVal()) < 11.5 ) {
	storePrefitResults(prefitParamMP2_[iMomBin], momMin, momMax, prefitResultMP2_);
	storePrefitResults(prefitParamWidth2_[iMomBin], momMin, momMax, prefitResultWidth2_);
      }
      storePrefitResults(prefitParamX0_[iMomBin], momMin, momMax, prefitResultX0_);
      storePrefitResults(prefitParamDeltaX1_[iMomBin], momMin, momMax, prefitResultDeltaX1_);

      canvas->Clear();
      canvas->SetLogy();
    
      TString frameTitle = Form("%s %s: P_{T} = %2.0f..%2.0f GeV", sep_->GetName(), decayMode_string.Data(), momMin, momMax);
      RooPlot* frame = sepTimesMom_prefit->frame(RooFit::Title(frameTitle.Data()), RooFit::Bins(120));

      datahist->plotOn(frame);
    
      prefitModel_[iMomBin]->plotOn(frame, RooFit::Range("combined"), RooFit::NormRange("combined"));

      frame->Draw();

      TLegend legend(0.70, 0.72, 0.89, 0.89, "", "brNDC");
      legend.SetBorderSize(0);
      legend.SetFillColor(0);
      TObject* graphModel = 0;
      TObject* graphData  = 0;
      for ( int iItem = 0; iItem < frame->numItems(); ++iItem ) {
	TString itemName = frame->nameOf(iItem);
	TObject* item = frame->findObject(itemName.Data());
	if      ( itemName.Contains("pdf")      ) graphModel = item;
	else if ( itemName.Contains("datahist") ) graphData  = item; 
      }
      if ( graphModel != 0 ) legend.AddEntry(graphModel, "Fit",    "l");
      if ( graphData  != 0 ) legend.AddEntry(graphData,  "TAUOLA", "p");
      legend.Draw();

      TString outputFileName_i = outputFileName;
      outputFileName_i = outputFileName_i.ReplaceAll(".", TString("_").Append(momBinName).Append("."));
      
      canvas->Update();
  
      canvas->SaveAs(outputFileName_i.Data());

      delete datahist;
    }

    delete canvas;

    delete inputFile;
   
    delete sepTimesMom_prefit;

    fitTF1(momParamGMean_, prefitResultGMean_, prefitCoeffValGMean_, prefitCoeffErrGMean_);
    fitTF1(momParamGSigma_, prefitResultGSigma_, prefitCoeffValGSigma_, prefitCoeffErrGSigma_);
    fitTF1(momParamAlpha_, prefitResultAlpha_, prefitCoeffValAlpha_, prefitCoeffErrAlpha_);
    fitTF1(momParamSlope_, prefitResultSlope_, prefitCoeffValSlope_, prefitCoeffErrSlope_);
    fitTF1(momParamOffset_, prefitResultOffset_, prefitCoeffValOffset_, prefitCoeffErrOffset_);
    fitTF1(momParamC_, prefitResultC_, prefitCoeffValC_, prefitCoeffErrC_);
    fitTF1(momParamMP1_, prefitResultMP1_, prefitCoeffValMP1_, prefitCoeffErrMP1_);
    fitTF1(momParamWidth1_, prefitResultWidth1_, prefitCoeffValWidth1_, prefitCoeffErrWidth1_);
    fitTF1(momParamMP2_, prefitResultMP2_, prefitCoeffValMP2_, prefitCoeffErrMP2_);
    fitTF1(momParamWidth2_, prefitResultWidth2_, prefitCoeffValWidth2_, prefitCoeffErrWidth2_);
    fitTF1(momParamX0_, prefitResultX0_, prefitCoeffValX0_, prefitCoeffErrX0_);
    fitTF1(momParamDeltaX1_, prefitResultDeltaX1_, prefitCoeffValDeltaX1_, prefitCoeffErrDeltaX1_);
  }

  RooAbsReal* buildMomDependentFitParameter(const TString& name, std::vector<RooRealVar*>& fitParameterCoefficients, 
					    std::vector<double>& prefitParamCoeffValues, std::vector<double>& prefitParamCoeffErrors,
					    const TString& momParametrizationFormula)
  {
    std::cout << "<fitManager::buildMomDependentFitParameter>:" << std::endl;

    TString decayMode_string = getDecayMode_string(decayMode_);

    TObjArray dependents;
  
    unsigned numCoefficients = prefitParamCoeffValues.size();
    std::cout << " numCoefficients = " << numCoefficients << std::endl;

    for ( unsigned iCoeff = 0; iCoeff < numCoefficients; ++iCoeff ) {
      TString fitParameterCoeffName = 
	Form("%s_%s_%s_%s_%s_p%u", name.Data(), mom_->GetName(), sep_->GetName(), decayMode_string.Data(), label_.Data(), iCoeff);
      std::cout << "--> creating fitParameterCoeff: name = " << fitParameterCoeffName.Data() << "," 
		<< " startValue = " << prefitParamCoeffValues[iCoeff] << " +/- " << prefitParamCoeffErrors[iCoeff] << std::endl;
      RooRealVar* fitParameterCoeff = 
	new RooRealVar(fitParameterCoeffName.Data(), fitParameterCoeffName.Data(), prefitParamCoeffValues[iCoeff]);
      dependents.Add(fitParameterCoeff);
      fitParameterCoefficients[iCoeff] = fitParameterCoeff;
    }

    dependents.Add(mom_);
  
    TString fitParameterName = 
      Form("%s_%s_%s_%s_%s", name.Data(), mom_->GetName(), sep_->GetName(), decayMode_string.Data(), label_.Data());
    TString fitParameterFormula = convertTFormulaToRooFitSyntax(momParametrizationFormula);
    std::cout << "--> creating momentum-dependent fitParameter: name = " << fitParameterName.Data() << ","
	      << " formula = " << fitParameterFormula.Data() << std::endl;
    RooAbsReal* fitParameter = new RooFormulaVar(fitParameterName.Data(), fitParameterFormula.Data(), RooArgList(dependents));
    
    return fitParameter;
  }

  void buildFitModel(const TString& outputFileName)
  {
    std::cout << "<fitManager::buildFitModel>:" << std::endl;

    printTimeStamp("<fitManager::buildFitModel>");

    TString decayMode_string = getDecayMode_string(decayMode_);
    
    fitParamGMean_   = buildMomDependentFitParameter("gmean", fitParamCoeffGMean_, 
						     prefitCoeffValGMean_, prefitCoeffErrGMean_, momParametrizationFormulaGMean_);
    fitParamGSigma_  = buildMomDependentFitParameter("gsigma", fitParamCoeffGSigma_, 
						     prefitCoeffValGSigma_, prefitCoeffErrGSigma_, momParametrizationFormulaGSigma_);
    //fitParamAlpha_   = buildMomDependentFitParameter("alpha", fitParamCoeffAlpha_, 
    //						       prefitCoeffValAlpha_, prefitCoeffErrAlpha_, momParametrizationFormulaAlpha_);
    fitParamAlpha_   = buildConstPrefitParameter("alpha", "AllMom", 0.);
    fitParamSlope_   = buildMomDependentFitParameter("slope", fitParamCoeffSlope_, 
						     prefitCoeffValSlope_, prefitCoeffErrSlope_, momParametrizationFormulaSlope_);
    //fitParamSlope_   = buildConstPrefitParameter("slope", "AllMom", 0.);
    //fitParamOffset_  = buildMomDependentFitParameter("offset", fitParamCoeffOffset_, 
    //						       prefitCoeffValOffset_, prefitCoeffErrOffset_, momParametrizationFormulaOffset_);
    fitParamOffset_  = buildConstPrefitParameter("offset", "AllMom", 0.);
    fitParamC_       = buildMomDependentFitParameter("C", fitParamCoeffC_, 
						     prefitCoeffValC_, prefitCoeffErrC_, momParametrizationFormulaC_);
    //fitParamC_       = buildConstPrefitParameter("C", "AllMom", 1.);
    fitParamMP1_     = buildMomDependentFitParameter("mp1", fitParamCoeffMP1_, 
						     prefitCoeffValMP1_, prefitCoeffErrMP1_, momParametrizationFormulaMP1_);
    fitParamWidth1_  = buildMomDependentFitParameter("width1", fitParamCoeffWidth1_, 
						     prefitCoeffValWidth1_, prefitCoeffErrWidth1_, momParametrizationFormulaWidth1_);
    fitParamMP2_     = buildMomDependentFitParameter("mp2", fitParamCoeffMP2_, 
						     prefitCoeffValMP2_, prefitCoeffErrMP2_, momParametrizationFormulaMP2_);
    fitParamWidth2_  = buildMomDependentFitParameter("width2", fitParamCoeffWidth2_,
						     prefitCoeffValWidth2_, prefitCoeffErrWidth2_, momParametrizationFormulaWidth2_);    
    fitParamX0_      = buildMomDependentFitParameter("x0", fitParamCoeffX0_, 
						     prefitCoeffValX0_, prefitCoeffErrX0_, momParametrizationFormulaX0_);
    fitParamDeltaX1_ = buildMomDependentFitParameter("dx1", fitParamCoeffDeltaX1_, 
						     prefitCoeffValDeltaX1_, prefitCoeffErrDeltaX1_, momParametrizationFormulaDeltaX1_);
    
    TString modelName = Form("pdf_%s_%s_%s_%s", decayMode_string.Data(), "AllMom", sep_->GetName(), label_.Data());
    fitModel_ = 
      new TauDecayKinePdf(modelName.Data(), modelName.Data(), 
			  *sepTimesMom_, 
			  *fitParamGMean_, *fitParamGSigma_, *fitParamAlpha_, *fitParamSlope_, *fitParamOffset_, *fitParamC_,
			  *fitParamMP1_, *fitParamWidth1_, *fitParamMP2_, *fitParamWidth2_, *fitParamX0_, *fitParamDeltaX1_);
  
    std::cout << "--> saving prefit results..." << std::endl;
    RooWorkspace* ws_prefit = new RooWorkspace("ws_prefit", "workspace");
    ws_prefit->import(*fitModel_);
    TString wsOutputFileName_prefit = outputFileName;
    wsOutputFileName_prefit.ReplaceAll("plots/fitTauDecayKinePlots", "mcTauDecayKine");
    wsOutputFileName_prefit.ReplaceAll(".eps", "_ws_prefit.root");
    std::cout << " wsOutputFileName_prefit = " << wsOutputFileName_prefit.Data() << std::endl;
    ws_prefit->writeToFile(wsOutputFileName_prefit.Data());
    std::cout << " done." << std::endl;

    std::cout << "--> estimate PDF timing (analytic integration enabled)..." << std::endl;
    pdfTimingTest(fitModel_, mom_, sepTimesMom_);
    std::cout << " done." << std::endl;
    
    fitModel_->disableAnalyticIntegration();
    
    std::cout << "--> estimate PDF timing (analytic integration disabled)..." << std::endl;
    pdfTimingTest(fitModel_, mom_, sepTimesMom_);
    std::cout << " done." << std::endl;
  }

  void runFit(RooAbsData* dataset, const TString& outputFileName)
  {
    std::cout << "<fitManager::runFit>:" << std::endl;

    printTimeStamp("<fitManager::runFit>");

    TString decayMode_string = getDecayMode_string(decayMode_);

    RooDataHist dataset_binned("dataset_binned", "dataset_binned", RooArgSet(*mom_, *sepTimesMom_), *dataset);
    
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

    fitModel_->enableAnalyticIntegration();

    std::cout << "--> saving fit results..." << std::endl;
    RooWorkspace* ws_fit = new RooWorkspace("ws_fit", "workspace");
    ws_fit->import(*fitModel_);
    TString wsOutputFileName_fit = outputFileName;
    wsOutputFileName_fit.ReplaceAll("plots/fitTauDecayKinePlots", "mcTauDecayKine");
    wsOutputFileName_fit.ReplaceAll(".eps", "_ws_fit.root");
    std::cout << " wsOutputFileName_fit = " << wsOutputFileName_fit.Data() << std::endl;
    ws_fit->writeToFile(wsOutputFileName_fit.Data());
    std::cout << " done." << std::endl;

    storeFitResults(fitParamCoeffGMean_, fitCoeffValGMean_, fitCoeffErrGMean_);
    storeFitResults(fitParamCoeffGSigma_, fitCoeffValGSigma_, fitCoeffErrGSigma_);
    storeFitResults(fitParamCoeffAlpha_, fitCoeffValAlpha_, fitCoeffErrAlpha_);
    storeFitResults(fitParamCoeffSlope_, fitCoeffValSlope_, fitCoeffErrSlope_);
    storeFitResults(fitParamCoeffOffset_, fitCoeffValOffset_, fitCoeffErrOffset_);
    storeFitResults(fitParamCoeffC_, fitCoeffValC_, fitCoeffErrC_);
    storeFitResults(fitParamCoeffMP1_, fitCoeffValMP1_, fitCoeffErrMP1_);
    storeFitResults(fitParamCoeffWidth1_, fitCoeffValWidth1_, fitCoeffErrWidth1_);
    storeFitResults(fitParamCoeffMP2_, fitCoeffValMP2_, fitCoeffErrMP2_);
    storeFitResults(fitParamCoeffWidth2_, fitCoeffValWidth2_, fitCoeffErrWidth2_);
    storeFitResults(fitParamCoeffX0_, fitCoeffValX0_, fitCoeffErrX0_);
    storeFitResults(fitParamCoeffDeltaX1_, fitCoeffValDeltaX1_, fitCoeffErrDeltaX1_);

    std::cout << "--> making control plots..." << std::endl;    
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

      TString frameTitle = Form("%s %s: P_{T} = %2.0f..%2.0f GeV", sep_->GetName(), decayMode_string.Data(), momMin, momMax);
      RooPlot* frame = sepTimesMom_->frame(RooFit::Title(frameTitle.Data()), RooFit::Bins(120));
      
      TString momSliceName = Form("momSlice%2.0fto%2.0f", momMin, momMax);
      std::stringstream momSliceCut;
      momSliceCut << momMin << " < " << mom_->GetName() << " && " << mom_->GetName() << " < " << momMax;
      std::cout << " momSliceCut = " << momSliceCut.str() << std::endl;
      
      RooAbsData* momSliceDataset = dataset->reduce(momSliceCut.str().data());
      std::cout << "--> Selected " << momSliceDataset->numEntries() << " TTree entries" 
		<< " in slice " << Form("mom = %2.0f..%2.0f GeV", momMin, momMax) << "..." << std::endl;
      RooDataHist momSliceDataset_binned("momSliceDataset_binned", "momSliceDataset_binned", 
					 RooArgSet(*mom_, *sepTimesMom_), *momSliceDataset);
      momSliceDataset_binned.plotOn(frame);
      //momSliceDataset->plotOn(frame);
     
      fitModel_->plotOn(frame, RooFit::ProjWData(RooArgSet(*mom_), momSliceDataset_binned));
      //fitModel_->plotOn(frame, RooFit::ProjWData(RooArgSet(*mom_), *momSliceDataset));

      frame->Draw();

      TLegend legend(0.70, 0.72, 0.89, 0.89, "", "brNDC");
      legend.SetBorderSize(0);
      legend.SetFillColor(0);
      TObject* graphModel = 0;
      TObject* graphData  = 0;
      for ( int iItem = 0; iItem < frame->numItems(); ++iItem ) {
	TString itemName = frame->nameOf(iItem);
	TObject* item = frame->findObject(itemName.Data());
	if      ( itemName.Contains("pdf")      ) graphModel = item;
	else if ( itemName.Contains("datahist") ) graphData  = item; 
      }
      if ( graphModel != 0 ) legend.AddEntry(graphModel, "Fit",    "l");
      if ( graphData  != 0 ) legend.AddEntry(graphData,  "TAUOLA", "p");
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

  RooRealVar* mom_;
  RooRealVar* sep_;
  RooRealVar* sepTimesMom_;
  int decayMode_;
  TString label_;
  TArrayD momBinning_;

  std::vector<RooRealVar*> prefitParamGMean_;
  std::vector<RooRealVar*> prefitParamGSigma_;
  std::vector<RooRealVar*> prefitParamAlpha_;
  std::vector<RooRealVar*> prefitParamSlope_;
  std::vector<RooRealVar*> prefitParamOffset_;
  std::vector<RooRealVar*> prefitParamC_;
  std::vector<RooRealVar*> prefitParamMP1_;
  std::vector<RooRealVar*> prefitParamWidth1_;
  std::vector<RooRealVar*> prefitParamMP2_;
  std::vector<RooRealVar*> prefitParamWidth2_;
  std::vector<RooRealVar*> prefitParamX0_;
  std::vector<RooRealVar*> prefitParamDeltaX1_;
  std::vector<TauDecayKinePdf*> prefitModel_;

  TGraphErrors* prefitResultGMean_;
  TGraphErrors* prefitResultGSigma_;
  TGraphErrors* prefitResultAlpha_;
  TGraphErrors* prefitResultSlope_;
  TGraphErrors* prefitResultOffset_;
  TGraphErrors* prefitResultC_;
  TGraphErrors* prefitResultMP1_;
  TGraphErrors* prefitResultWidth1_;
  TGraphErrors* prefitResultMP2_;
  TGraphErrors* prefitResultWidth2_;
  TGraphErrors* prefitResultX0_;
  TGraphErrors* prefitResultDeltaX1_;

  TString momParametrizationFormulaGMean_;
  TString momParametrizationFormulaGSigma_;
  TString momParametrizationFormulaAlpha_;
  TString momParametrizationFormulaSlope_;
  TString momParametrizationFormulaOffset_;
  TString momParametrizationFormulaC_;
  TString momParametrizationFormulaMP1_;
  TString momParametrizationFormulaWidth1_;
  TString momParametrizationFormulaMP2_;
  TString momParametrizationFormulaWidth2_;
  TString momParametrizationFormulaX0_;
  TString momParametrizationFormulaDeltaX1_;

  TF1* momParamGMean_;
  TF1* momParamGSigma_;
  TF1* momParamAlpha_;
  TF1* momParamSlope_;
  TF1* momParamOffset_;
  TF1* momParamC_;
  TF1* momParamMP1_;
  TF1* momParamWidth1_;
  TF1* momParamMP2_;
  TF1* momParamWidth2_;
  TF1* momParamX0_;
  TF1* momParamDeltaX1_;

  std::vector<double> prefitCoeffValGMean_;
  std::vector<double> prefitCoeffErrGMean_;
  std::vector<double> prefitCoeffValGSigma_;
  std::vector<double> prefitCoeffErrGSigma_;
  std::vector<double> prefitCoeffValAlpha_;
  std::vector<double> prefitCoeffErrAlpha_;
  std::vector<double> prefitCoeffValSlope_;
  std::vector<double> prefitCoeffErrSlope_;
  std::vector<double> prefitCoeffValOffset_;
  std::vector<double> prefitCoeffErrOffset_;
  std::vector<double> prefitCoeffValC_;
  std::vector<double> prefitCoeffErrC_;
  std::vector<double> prefitCoeffValMP1_;
  std::vector<double> prefitCoeffErrMP1_;
  std::vector<double> prefitCoeffValWidth1_;
  std::vector<double> prefitCoeffErrWidth1_;
  std::vector<double> prefitCoeffValMP2_;
  std::vector<double> prefitCoeffErrMP2_;
  std::vector<double> prefitCoeffValWidth2_;
  std::vector<double> prefitCoeffErrWidth2_;
  std::vector<double> prefitCoeffValX0_;
  std::vector<double> prefitCoeffErrX0_;
  std::vector<double> prefitCoeffValDeltaX1_;
  std::vector<double> prefitCoeffErrDeltaX1_;
  
  std::vector<RooRealVar*> fitParamCoeffGMean_;
  RooAbsReal* fitParamGMean_;
  std::vector<RooRealVar*> fitParamCoeffGSigma_;
  RooAbsReal* fitParamGSigma_;
  std::vector<RooRealVar*> fitParamCoeffAlpha_;
  RooAbsReal* fitParamAlpha_;
  std::vector<RooRealVar*> fitParamCoeffSlope_;
  RooAbsReal* fitParamSlope_;
  std::vector<RooRealVar*> fitParamCoeffOffset_;
  RooAbsReal* fitParamOffset_;
  std::vector<RooRealVar*> fitParamCoeffC_;
  RooAbsReal* fitParamC_;
  std::vector<RooRealVar*> fitParamCoeffMP1_;
  RooAbsReal* fitParamMP1_;
  std::vector<RooRealVar*> fitParamCoeffWidth1_;
  RooAbsReal* fitParamWidth1_;
  std::vector<RooRealVar*> fitParamCoeffMP2_;
  RooAbsReal* fitParamMP2_;
  std::vector<RooRealVar*> fitParamCoeffWidth2_;
  RooAbsReal* fitParamWidth2_;
  std::vector<RooRealVar*> fitParamCoeffX0_;
  RooAbsReal* fitParamX0_;
  std::vector<RooRealVar*> fitParamCoeffDeltaX1_;
  RooAbsReal* fitParamDeltaX1_;
  TauDecayKinePdf* fitModel_;

  std::vector<double> fitCoeffValGMean_;
  std::vector<double> fitCoeffErrGMean_;
  std::vector<double> fitCoeffValGSigma_;
  std::vector<double> fitCoeffErrGSigma_;
  std::vector<double> fitCoeffValAlpha_;
  std::vector<double> fitCoeffErrAlpha_;
  std::vector<double> fitCoeffValSlope_;
  std::vector<double> fitCoeffErrSlope_;
  std::vector<double> fitCoeffValOffset_;
  std::vector<double> fitCoeffErrOffset_;
  std::vector<double> fitCoeffValC_;
  std::vector<double> fitCoeffErrC_;
  std::vector<double> fitCoeffValMP1_;
  std::vector<double> fitCoeffErrMP1_;
  std::vector<double> fitCoeffValWidth1_;
  std::vector<double> fitCoeffErrWidth1_;
  std::vector<double> fitCoeffValMP2_;
  std::vector<double> fitCoeffErrMP2_;
  std::vector<double> fitCoeffValWidth2_;
  std::vector<double> fitCoeffErrWidth2_;
  std::vector<double> fitCoeffValX0_;
  std::vector<double> fitCoeffErrX0_;
  std::vector<double> fitCoeffValDeltaX1_;
  std::vector<double> fitCoeffErrDeltaX1_;
};

//
//-------------------------------------------------------------------------------
//

int main(int argc, const char* argv[])
{
  if ( argc < 5 ) {
    std::cerr << "Usage: ./fitTauDecayKine inputFileNames decayMode leg parametrization selection" << std::endl;
    return 1;
  }

  printTimeStamp("<fitTauDecayKine::main (begin)>");
  
  TString inputFileNames_ntuple = argv[1];
  //inputFileNames_ntuple = "/data2/friis/PtBalanceNtupleData_v4/ptBalanceData_*_ggAH.root";
  //inputFileNames_ntuple = "/data2/friis/PtBalanceNtupleData_v4/ptBalanceData_mass_200_ggAH.root";

  std::vector<unsigned> decayModesToRun;
  for ( unsigned iDecayMode = kElectron_Muon; iDecayMode <= kThreeProng1Pi0; ++iDecayMode ) {
    if ( getDecayMode_string(iDecayMode) == argv[2] ) decayModesToRun.push_back(iDecayMode);
  }

  if ( decayModesToRun.size() == 0 ) 
    throw cms::Exception("fitTauDecayKine")
      << "Invalid Configuration Parameter 'decayMode' = " << argv[2] << " !!\n";

  bool runLeg1     = ( std::string(argv[3]) == "leg1"     ) ? true : false;
  bool runLeg2     = ( std::string(argv[3]) == "leg2"     ) ? true : false;

  bool runDeltaR   = ( std::string(argv[4]) == "dR"       ) ? true : false;
  bool runAngle    = ( std::string(argv[4]) == "angle"    ) ? true : false;
  if ( !(runDeltaR || runAngle) )
    throw cms::Exception("fitTauDecayKine")
      << "Invalid Configuration Parameter 'parametrization' = " << argv[4] << " !!\n";

  bool runAll      = ( std::string(argv[5]) == "all"      ) ? true : false;
  bool runSelected = ( std::string(argv[5]) == "selected" ) ? true : false;
  if ( !(runAll || runSelected) )
    throw cms::Exception("fitTauDecayKine")
      << "Invalid Configuration Parameter 'selection' = " << argv[5] << " !!\n";

  bool runFit = false;
  //bool runFit = true;
  
  TString inputFileName_histograms = "makeTauDecayKinePlots.root";

  double decayModeElectron_Muon_encoded = encodeDecayMode(kElectron_Muon);
  double decayModeHadMin_encoded = encodeDecayMode(kOneProng0Pi0);
  double decayModeHadMax_encoded = encodeDecayMode(kThreeProng1Pi0);

  TString nanFilter;
  nanFilter.Append("!(TMath::IsNaN(leg1Pt) || TMath::IsNaN(leg1VisPt) ||");
  nanFilter.Append("  TMath::IsNaN(leg1VisInvisAngleLab) || TMath::IsNaN(leg1VisInvisDeltaRLab) ||");
  nanFilter.Append("  TMath::IsNaN(leg2Pt) || TMath::IsNaN(leg2VisPt) ||");
  nanFilter.Append("  TMath::IsNaN(leg2VisInvisAngleLab) || TMath::IsNaN(leg2VisInvisDeltaRLab))");

  TString visMomCuts_leg1;
  visMomCuts_leg1.Append(Form("(leg1VisPt > 15. && TMath::Abs(leg1VisEta) < 2.1 &&"));
  visMomCuts_leg1.Append(Form(" leg1DecayMode > %2.1f && leg1DecayMode < %2.1f) || ",
			      decayModeElectron_Muon_encoded - 0.1, decayModeElectron_Muon_encoded + 0.1));
  visMomCuts_leg1.Append(Form("(leg1VisPt > 20. && TMath::Abs(leg1VisEta) < 2.3 &&"));
  visMomCuts_leg1.Append(Form(" leg1DecayMode > %2.1f && leg1DecayMode < %2.1f)", 
			      decayModeHadMin_encoded - 0.1, decayModeHadMax_encoded + 0.1));
  std::cout << "visMomCuts_leg1 = " << visMomCuts_leg1.Data() << std::endl;
  TString visMomCuts_leg2;
  visMomCuts_leg2.Append(Form("(leg2VisPt > 15. && TMath::Abs(leg2VisEta) < 2.1 &&"));
  visMomCuts_leg2.Append(Form(" leg2DecayMode > %2.1f && leg2DecayMode < %2.1f) || ",
			      decayModeElectron_Muon_encoded - 0.1, decayModeElectron_Muon_encoded + 0.1));
  visMomCuts_leg2.Append(Form("(leg2VisPt > 20. && TMath::Abs(leg2VisEta) < 2.3 &&"));
  visMomCuts_leg2.Append(Form(" leg2DecayMode > %2.1f && leg2DecayMode < %2.1f)", 
			      decayModeHadMin_encoded - 0.1, decayModeHadMax_encoded + 0.1));
  std::cout << "visMomCuts_leg2 = " << visMomCuts_leg2.Data() << std::endl;
  
//--- stop RooT from keeping references to all histograms
  TH1::AddDirectory(false); 

  gROOT->SetBatch(true);

  // Load the data
  std::cout << "Loading data from " << inputFileNames_ntuple <<  std::endl;
  TChain* dataTree = 0;
  if ( runFit ) {
    dataTree = new TChain("makePtBalanceNtuple/ptBalanceNtuple");
    dataTree->Add(inputFileNames_ntuple);
  }

  /*****************************************************************************
   ********                   Define variabes and data                      ****
   *****************************************************************************/

  RooRealVar leg1Energy("leg1Energy", "leg1Energy", 0., 350.);
  RooRealVar leg1Pt("leg1Pt", "leg1Pt", 0., 250.);
  RooRealVar leg1VisPt("leg1VisPt", "leg1VisPt", 0., 200.);
  RooRealVar leg1VisEta("leg1VisEta", "leg1VisEta", -2.5, +2.5);
  RooRealVar leg1DecayMode("leg1DecayMode", "leg1DecayMode", -1.5, +11.5);
  RooRealVar leg1VisInvisAngleLab("leg1VisInvisAngleLab", "leg1VisInvisAngleLab", 0., 1.);
  RooRealVar leg1VisInvisDeltaRLab("leg1VisInvisDeltaRLab", "leg1VisInvisDeltaRLab", 0., 1.);
  RooRealVar leg1VisInvisAngleLabTimesEnergy("leg1VisInvisAngleLabTimesEnergy", "leg1VisInvisAngleLabTimesEnergy", 0., 12.);
  RooRealVar leg1VisInvisDeltaRLabTimesPt("leg1VisInvisDeltaRLabTimesPt", "leg1VisInvisDeltaRLabTimesPt", 0., 12.);

  leg1Energy.setBins(350);
  leg1Pt.setBins(250);
  leg1VisInvisAngleLab.setBins(120);
  leg1VisInvisDeltaRLab.setBins(120);
  leg1VisInvisAngleLabTimesEnergy.setBins(120);
  leg1VisInvisDeltaRLabTimesPt.setBins(120);

  RooRealVar leg2Energy("leg2Energy", "leg2Energy", 0., 350.);
  RooRealVar leg2Pt("leg2Pt", "leg2Pt", 0., 250.);
  RooRealVar leg2VisPt("leg2VisPt", "leg2VisPt", 0., 200.);
  RooRealVar leg2VisEta("leg2VisEta", "leg2VisEta", -2.5, +2.5);
  RooRealVar leg2DecayMode("leg2DecayMode", "leg2DecayMode", -1.5, +11.5);
  RooRealVar leg2VisInvisAngleLab("leg2VisInvisAngleLab", "leg2VisInvisAngleLab", 0., 1.);
  RooRealVar leg2VisInvisDeltaRLab("leg2VisInvisDeltaRLab", "leg2VisInvisDeltaRLab", 0., 1.);
  RooRealVar leg2VisInvisAngleLabTimesEnergy("leg2VisInvisAngleLabTimesEnergy", "leg2VisInvisAngleLabTimesEnergy", 0., 12.);
  RooRealVar leg2VisInvisDeltaRLabTimesPt("leg2VisInvisDeltaRLabTimesPt", "leg2VisInvisDeltaRLabTimesPt", 0., 12.);  

  leg2Energy.setBins(350);
  leg2Pt.setBins(250);
  leg2VisInvisAngleLab.setBins(120);
  leg2VisInvisDeltaRLab.setBins(120);
  leg2VisInvisAngleLabTimesEnergy.setBins(120);
  leg2VisInvisDeltaRLabTimesPt.setBins(120);

  TObjArray variables;
  variables.Add(&leg1Energy);
  variables.Add(&leg1Pt);
  variables.Add(&leg1VisPt);
  variables.Add(&leg1VisEta);
  variables.Add(&leg1DecayMode);
  variables.Add(&leg1VisInvisAngleLab);
  variables.Add(&leg1VisInvisDeltaRLab);
  variables.Add(&leg1VisInvisAngleLabTimesEnergy);
  variables.Add(&leg1VisInvisDeltaRLabTimesPt);
  variables.Add(&leg2Energy);
  variables.Add(&leg2Pt);
  variables.Add(&leg2VisPt);
  variables.Add(&leg2VisEta);
  variables.Add(&leg2DecayMode);
  variables.Add(&leg2VisInvisAngleLab);
  variables.Add(&leg2VisInvisDeltaRLab);
  variables.Add(&leg2VisInvisAngleLabTimesEnergy);
  variables.Add(&leg2VisInvisDeltaRLabTimesPt);

  RooDataSet* dataset_all           = 0;
  RooAbsData* dataset_notNaN        = 0;
  RooAbsData* dataset_selected_leg1 = 0;
  RooAbsData* dataset_selected_leg2 = 0;
  if ( runFit ) {
    dataset_all = new RooDataSet("dataset", "datasetset", RooArgSet(variables), RooFit::Import(*dataTree));
    std::cout << "Processing " << dataset_all->numEntries() << " TTree entries..." << std::endl;
    
    dataset_notNaN = dataset_all->reduce(nanFilter.Data());
    std::cout << "--> Passing anti-NaN filter = " << dataset_notNaN->numEntries() << std::endl;

    dataset_selected_leg1 = dataset_notNaN->reduce(visMomCuts_leg1.Data());
    std::cout << "--> Passing vis. Momentum Cuts for leg1 = " << dataset_selected_leg1->numEntries() << std::endl;
    dataset_selected_leg2 = dataset_notNaN->reduce(visMomCuts_leg2.Data());
    std::cout << "--> Passing vis. Momentum Cuts for leg2 = " << dataset_selected_leg2->numEntries() << std::endl;
  }

  std::map<std::string, fitManager*> fitResults;
  
  for ( std::vector<unsigned>::const_iterator decayModeToRun = decayModesToRun.begin();
	decayModeToRun != decayModesToRun.end(); ++decayModeToRun ) {
    double decayMode_encoded = encodeDecayMode(*decayModeToRun);
  
    TString decayMode_string = getDecayMode_string(*decayModeToRun);
  
    TString leg1DecayModeSelection = Form("leg1DecayMode > %2.1f && leg1DecayMode < %2.1f",
					  decayMode_encoded - 0.1, decayMode_encoded + 0.1);
    if ( (*decayModeToRun) == kOneProngGt0Pi0 ) 
      leg1DecayModeSelection = 
	Form("leg1DecayMode > %2.1f && leg1DecayMode < %2.1f", 
	     encodeDecayMode(kOneProng1Pi0) - 0.1, encodeDecayMode(kOneProng2Pi0) + 0.1);

    TString leg2DecayModeSelection = Form("leg2DecayMode > %2.1f && leg2DecayMode < %2.1f",
					  decayMode_encoded - 0.1, decayMode_encoded + 0.1);
    if ( (*decayModeToRun) == kOneProngGt0Pi0 ) 
      leg2DecayModeSelection = 
	Form("leg2DecayMode > %2.1f && leg2DecayMode < %2.1f", 
	     encodeDecayMode(kOneProng1Pi0) - 0.1, encodeDecayMode(kOneProng2Pi0) + 0.1);
  
    TString selection_string = "";
    if ( runAll      ) selection_string.Append("_all");
    if ( runSelected ) selection_string.Append("_selected");
      
    TArrayD momBinning = getBinningMom(decayMode_string, selection_string);
    TArrayD sepBinning = getBinningSepTimesMom(decayMode_string, selection_string);
    
    /*****************************************************************************
     ******** Run fit for leg1 with no cuts on visible decay products applied ****
     *****************************************************************************/

    if ( runLeg1 && runAll && runDeltaR ) {
      TString fitManagerName_leg1_dR_all = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg1", "dR", "all");
      if ( fitResults.find(fitManagerName_leg1_dR_all.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg1_dR_all.Data()] = 
	  new fitManager(&leg1Pt, &leg1VisInvisDeltaRLab, &leg1VisInvisDeltaRLabTimesPt, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg1", "dR", "all"), momBinning);
      }
      fitManager* fitManager_leg1_dR_all = fitResults[fitManagerName_leg1_dR_all.Data()];
      TString outputFileName_leg1_dR_all = Form("plots/fitTauDecayKinePlots_%s_leg1_dR_all.eps", decayMode_string.Data());
      fitManager_leg1_dR_all->runPrefit(inputFileName_histograms, 
					"VisInvisDeltaRLabTimesPt", "all", outputFileName_leg1_dR_all);      
      fitManager_leg1_dR_all->buildFitModel(outputFileName_leg1_dR_all);
      if ( runFit ) { 
	RooAbsData* dataset_leg1_all = dataset_notNaN->reduce(leg1DecayModeSelection.Data());
	fitManager_leg1_dR_all->runFit(dataset_leg1_all, outputFileName_leg1_dR_all);
      }
    }

    if ( runLeg1 && runAll && runAngle ) {
      TString fitManagerName_leg1_angle_all = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg1", "angle", "all");
      if ( fitResults.find(fitManagerName_leg1_angle_all.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg1_angle_all.Data()] = 
	  new fitManager(&leg1Energy, &leg1VisInvisAngleLab, &leg1VisInvisAngleLabTimesEnergy, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg1", "angle", "all"), momBinning);
      }
      fitManager* fitManager_leg1_angle_all = fitResults[fitManagerName_leg1_angle_all.Data()];
      TString outputFileName_leg1_angle_all = 
	Form("plots/fitTauDecayKinePlots_%s_leg1_angle_all.eps", decayMode_string.Data());
      fitManager_leg1_angle_all->runPrefit(inputFileName_histograms, 
					   "VisInvisAngleLabTimesEnergy", "all", outputFileName_leg1_angle_all);
      fitManager_leg1_angle_all->buildFitModel(outputFileName_leg1_angle_all);
      if ( runFit ) {
	RooAbsData* dataset_leg1_all = dataset_notNaN->reduce(leg1DecayModeSelection.Data());
	fitManager_leg1_angle_all->runFit(dataset_leg1_all, outputFileName_leg1_angle_all);
      }
    }

    /*****************************************************************************
     ********   Run fit for leg1 with cuts on visible decay products applied  ****
     *****************************************************************************/

    if ( runLeg1 && runSelected && runDeltaR ) {
      TString fitManagerName_leg1_dR_selected = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg1", "dR", "selected");
      if ( fitResults.find(fitManagerName_leg1_dR_selected.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg1_dR_selected.Data()] = 
	  new fitManager(&leg1Pt, &leg1VisInvisDeltaRLab, &leg1VisInvisDeltaRLabTimesPt, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg1", "dR", "selected1"), momBinning);
      }
      fitManager* fitManager_leg1_dR_selected = fitResults[fitManagerName_leg1_dR_selected.Data()];
      TString outputFileName_leg1_dR_selected = Form("plots/fitTauDecayKinePlots_%s_leg1_dR_selected.eps", decayMode_string.Data());
      fitManager_leg1_dR_selected->runPrefit(inputFileName_histograms, 
					     "VisInvisDeltaRLabTimesPt", "selected1", outputFileName_leg1_dR_selected);
      fitManager_leg1_dR_selected->buildFitModel(outputFileName_leg1_dR_selected);
      if ( runFit ) {
	RooAbsData* dataset_leg1_selected = dataset_selected_leg1->reduce(leg1DecayModeSelection.Data());
	fitManager_leg1_dR_selected->runFit(dataset_leg1_selected, outputFileName_leg1_dR_selected);
      }
    }

    if ( runLeg1 && runSelected && runAngle ) {
      TString fitManagerName_leg1_angle_selected = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg1", "angle", "selected");
      if ( fitResults.find(fitManagerName_leg1_angle_selected.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg1_angle_selected.Data()] = 
	  new fitManager(&leg1Energy, &leg1VisInvisAngleLab, &leg1VisInvisAngleLabTimesEnergy, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg1", "angle", "selected1"), momBinning);
      }
      fitManager* fitManager_leg1_angle_selected = fitResults[fitManagerName_leg1_angle_selected.Data()];
      TString outputFileName_leg1_angle_selected = 
	Form("plots/fitTauDecayKinePlots_%s_leg1_angle_selected.eps", decayMode_string.Data());
      fitManager_leg1_angle_selected->runPrefit(inputFileName_histograms, 
						"VisInvisAngleLabTimesEnergy", "selected1", outputFileName_leg1_angle_selected);
      fitManager_leg1_angle_selected->buildFitModel(outputFileName_leg1_angle_selected);
      if ( runFit ) {
	RooAbsData* dataset_leg1_selected = dataset_selected_leg1->reduce(leg1DecayModeSelection.Data());
	fitManager_leg1_angle_selected->runFit(dataset_leg1_selected, outputFileName_leg1_angle_selected);
      }
    }

    /*****************************************************************************
     ******** Run fit for leg2 with no cuts on visible decay products applied ****
     *****************************************************************************/

    if ( runLeg2 && runAll && runDeltaR ) {
      TString fitManagerName_leg2_dR_all = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg2", "dR", "all");
      if ( fitResults.find(fitManagerName_leg2_dR_all.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg2_dR_all.Data()] = 
	  new fitManager(&leg2Pt, &leg2VisInvisDeltaRLab, &leg2VisInvisDeltaRLabTimesPt, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg2", "dR", "all"), momBinning);
      }
      fitManager* fitManager_leg2_dR_all = fitResults[fitManagerName_leg2_dR_all.Data()];
      TString outputFileName_leg2_dR_all = Form("plots/fitTauDecayKinePlots_%s_leg2_dR_all.eps", decayMode_string.Data());
      fitManager_leg2_dR_all->runPrefit(inputFileName_histograms, 
					"VisInvisDeltaRLabTimesPt", "all", outputFileName_leg2_dR_all);
      fitManager_leg2_dR_all->buildFitModel(outputFileName_leg2_dR_all);
      if ( runFit ) {
	RooAbsData* dataset_leg2_all = dataset_notNaN->reduce(leg2DecayModeSelection.Data());
	fitManager_leg2_dR_all->runFit(dataset_leg2_all, outputFileName_leg2_dR_all);
      }
    }

    if ( runLeg2 && runAll && runAngle ) {    
      TString fitManagerName_leg2_angle_all = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg2", "angle", "all");
      if ( fitResults.find(fitManagerName_leg2_angle_all.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg2_angle_all.Data()] = 
	  new fitManager(&leg2Energy, &leg2VisInvisAngleLab, &leg2VisInvisAngleLabTimesEnergy, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg2", "angle", "all"), momBinning);
      }
      fitManager* fitManager_leg2_angle_all = fitResults[fitManagerName_leg2_angle_all.Data()];
      TString outputFileName_leg2_angle_all = 
	Form("plots/fitTauDecayKinePlots_%s_leg2_angle_all.eps", decayMode_string.Data());
      fitManager_leg2_angle_all->runPrefit(inputFileName_histograms, 
					   "VisInvisAngleLabTimesEnergy", "all", outputFileName_leg2_angle_all);
      fitManager_leg2_angle_all->buildFitModel(outputFileName_leg2_angle_all);
      if ( runFit ) {
	RooAbsData* dataset_leg2_all = dataset_notNaN->reduce(leg2DecayModeSelection.Data());
	fitManager_leg2_angle_all->runFit(dataset_leg2_all, outputFileName_leg2_angle_all);
      }
    }
    
    /*****************************************************************************
     ********   Run fit for leg2 with cuts on visible decay products applied  ****
     *****************************************************************************/

    if ( runLeg2 && runSelected && runDeltaR ) {
      TString fitManagerName_leg2_dR_selected = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg2", "dR", "selected");
      if ( fitResults.find(fitManagerName_leg2_dR_selected.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg2_dR_selected.Data()] = 
	  new fitManager(&leg2Pt, &leg2VisInvisDeltaRLab, &leg2VisInvisDeltaRLabTimesPt, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg2", "dR", "selected1"), momBinning);
      }
      fitManager* fitManager_leg2_dR_selected = fitResults[fitManagerName_leg2_dR_selected.Data()];
      TString outputFileName_leg2_dR_selected = Form("plots/fitTauDecayKinePlots_%s_leg2_dR_selected.eps", decayMode_string.Data());
      fitManager_leg2_dR_selected->runPrefit(inputFileName_histograms, 
					     "VisInvisDeltaRLabTimesPt", "selected1", outputFileName_leg2_dR_selected);
      fitManager_leg2_dR_selected->buildFitModel(outputFileName_leg2_dR_selected);
      if ( runFit ) {
	RooAbsData* dataset_leg2_selected = dataset_selected_leg2->reduce(leg2DecayModeSelection.Data());
	fitManager_leg2_dR_selected->runFit(dataset_leg2_selected, outputFileName_leg2_dR_selected);
      }
    }
     
    if ( runLeg2 && runSelected && runAngle ) {
      TString fitManagerName_leg2_angle_selected = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg2", "angle", "selected");
      if ( fitResults.find(fitManagerName_leg2_angle_selected.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg2_angle_selected.Data()] = 
	  new fitManager(&leg2Energy, &leg2VisInvisAngleLab, &leg2VisInvisAngleLabTimesEnergy, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg2", "angle", "selected1"), momBinning);
      }
      fitManager* fitManager_leg2_angle_selected = fitResults[fitManagerName_leg2_angle_selected.Data()];
      TString outputFileName_leg2_angle_selected = 
	Form("plots/fitTauDecayKinePlots_%s_leg2_angle_selected.eps", decayMode_string.Data());
      fitManager_leg2_angle_selected->runPrefit(inputFileName_histograms, 
						"VisInvisAngleLabTimesEnergy", "selected1", outputFileName_leg2_angle_selected);
      fitManager_leg2_angle_selected->buildFitModel(outputFileName_leg2_angle_selected);
      if ( runFit ) {
	RooAbsData* dataset_leg2_selected = dataset_selected_leg2->reduce(leg2DecayModeSelection.Data());
	fitManager_leg2_angle_selected->runFit(dataset_leg2_selected, outputFileName_leg2_angle_selected);
      }
    }
  }

  delete dataTree;

  printTimeStamp("<fitTauDecayKine::main (end)>");
}
