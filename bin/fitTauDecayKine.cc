
#include "TauAnalysis/GenSimTools/bin/tauDecayKineAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/TauDecayKinePdf2.h"

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
    double sepTimesMomValue = rnd.Uniform( 0.,  25.);

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

double compFallingEdgePos(TH1* histogram, int decayMode, double& errEdgePosRight, double& errEdgePosLeft)
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

  if ( decayMode == kElectron_Muon  ||
       decayMode == kOneProngGt0Pi0 ||
       decayMode == kThreeProng0Pi0 ) {
    binFallingEdge -= 2; // CV: "phenomenological" correction...
  }

  if ( binFallingEdge > 0 ) {
    errEdgePosRight = TMath::Abs(histogram->GetBinContent(binFallingEdge)/minWindow5derrivative);
    std::cout << " errEdgePosRight = " << errEdgePosRight << std::endl;
    
    double histogramMax_y = histogram->GetBinContent(histogram->GetMaximumBin());      
    errEdgePosLeft  = TMath::Abs((histogramMax_y - histogram->GetBinContent(binFallingEdge))/minWindow5derrivative);
    std::cout << " errEdgePosLeft = " << errEdgePosLeft << std::endl;
  }

  return fallindEdgePos;
}

void compPrePeakPos(TH1* histogram, int decayMode, double& prePeakPosLeft, double& prePeakPosRight)
{
  int binMax = histogram->GetMaximumBin();
  double binContentMax = histogram->GetBinContent(binMax);

  bool hasMinPassed = false;

  int binPrePeakMax = -1;
  double binContentPrePeakMax = -1.;

  int binPrePeakMin = -1;
  double binContentPrePeakMin = binContentMax;

  for ( int iBin = binMax; iBin >= 1; --iBin ) {
    double binContent = histogram->GetBinContent(iBin);
    
    if ( !hasMinPassed && binContent < (0.5*binContentMax) ) {
      hasMinPassed = true;
      continue;
    }
    
    if (  hasMinPassed && binContent >  binContentPrePeakMax  ) {
      binPrePeakMax = iBin;
      binContentPrePeakMax = binContent;
    }
  }

  if ( binPrePeakMax != -1 ) {
    for ( int iBin = binMax; iBin >= binPrePeakMax; --iBin ) {
      double binContent = histogram->GetBinContent(iBin);
      
      if ( binContent < binContentPrePeakMin  ) {
	binPrePeakMin = iBin;
	binContentPrePeakMin = binContent;
      }
    }
  } else {
    binPrePeakMax = binMax;
    binPrePeakMin = binMax;
  }

  prePeakPosLeft  = histogram->GetBinCenter(binPrePeakMax);
  prePeakPosRight = histogram->GetBinCenter(binPrePeakMin);
}

double compFCN(TH1* histogram, int iBin0, int iBin1, double mp, double width)
{
  //std::cout << "<compFCN>:" << std::endl;

  double retVal = 0.;
  
  double histogram_sum = 0.;
  double landau_integral = 0.;
  for ( int iBin = iBin0; iBin <= iBin1; ++iBin ) {
    histogram_sum += histogram->GetBinContent(iBin);
    landau_integral += ::ROOT::Math::landau_pdf((histogram->GetBinCenter(iBin) - mp)/width);
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
    
    double landau = landau_norm*::ROOT::Math::landau_pdf((binCenter - mp)/width);
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
      prefitParamMP1_(momBinning.GetSize() - 1),
      prefitParamWidth1_(momBinning.GetSize() - 1),
      prefitParamGMean_(momBinning.GetSize() - 1),
      prefitParamGSigma_(momBinning.GetSize() - 1),
      prefitParamAlpha_(momBinning.GetSize() - 1),
      prefitParamSlope_(momBinning.GetSize() - 1),
      prefitParamOffset_(momBinning.GetSize() - 1),
      prefitParamC_(momBinning.GetSize() - 1),
      prefitParamMP2_(momBinning.GetSize() - 1),
      prefitParamWidth2_(momBinning.GetSize() - 1),
      prefitParamMP3_(momBinning.GetSize() - 1),
      prefitParamWidth3_(momBinning.GetSize() - 1),
      prefitParamX0_(momBinning.GetSize() - 1),
      prefitParamDeltaX1_(momBinning.GetSize() - 1),
      prefitParamDeltaX2_(momBinning.GetSize() - 1),
      prefitModel_(momBinning.GetSize() - 1),
      momParamMP1_(0),
      momParamWidth1_(0),
      momParamGMean_(0),
      momParamGSigma_(0),
      momParamAlpha_(0),
      momParamSlope_(0),
      momParamOffset_(0),
      momParamC_(0),
      momParamMP2_(0),
      momParamWidth2_(0),
      momParamMP3_(0),
      momParamWidth3_(0),
      momParamX0_(0),
      momParamDeltaX1_(0),
      momParamDeltaX2_(0),
      fitParamMP1_(0),
      fitParamWidth1_(0),
      fitParamGMean_(0),
      fitParamGSigma_(0),
      fitParamAlpha_(0),
      fitParamSlope_(0),
      fitParamOffset_(0),
      fitParamC_(0),
      fitParamMP2_(0),
      fitParamWidth2_(0),
      fitParamMP3_(0),
      fitParamWidth3_(0),
      fitParamX0_(0),
      fitParamDeltaX1_(0),
      fitParamDeltaX2_(0),
      fitModel_(0)
  {
    std::cout << "<fitManager::fitManager>:" << std::endl;
    std::cout << " mom = " << mom_->GetName() << std::endl;
    std::cout << " dR = " << sep_->GetName() << std::endl;
    std::cout << " decayMode = " << getDecayMode_string(decayMode_) << std::endl;

    isAngle_    = false;
    isDeltaR_   = false;
    if      ( label.Contains("_angle")    ) isAngle_    = true;
    else if ( label.Contains("_dR")       ) isDeltaR_   = true;
    else assert (0);

    isAll_      = false;
    isSelected_ = false;
    if      ( label.Contains("_all")      ) isAll_      = true;
    else if ( label.Contains("_selected") ) isSelected_ = true;
    else assert (0);

    prefitResultMP1_ = new TGraphErrors();
    prefitResultWidth1_ = new TGraphErrors();
    prefitResultGMean_ = new TGraphErrors();
    prefitResultGSigma_ = new TGraphErrors();
    prefitResultAlpha_ = new TGraphErrors();
    prefitResultSlope_ = new TGraphErrors();
    prefitResultOffset_ = new TGraphErrors();
    prefitResultC_ = new TGraphErrors();
    prefitResultMP2_ = new TGraphErrors();
    prefitResultWidth2_ = new TGraphErrors();
    prefitResultMP3_ = new TGraphErrors();
    prefitResultWidth3_ = new TGraphErrors();
    prefitResultX0_ = new TGraphErrors();
    prefitResultDeltaX1_ = new TGraphErrors();
    prefitResultDeltaX2_ = new TGraphErrors();
    
    momParametrizationFormulaMP1_ = "[0] + [1]*x";
    setMinMP1_ = false;
    setMaxMP1_ = false;
    fitParamCoeffMP1_.resize(2);
    momParamMP1_ = bookTF1("MP1", momParametrizationFormulaMP1_);
    momParametrizationFormulaWidth1_ = "[0] + [1]*x";
    setMinWidth1_ = false;
    setMaxWidth1_ = false;
    fitParamCoeffWidth1_.resize(2);
    momParamWidth1_ = bookTF1("Width1", momParametrizationFormulaWidth1_);
    if ( decayMode_ == kElectron_Muon && isAll_ ) {
      momParametrizationFormulaGMean_ = "[0] + [1]/TMath::Power(x - 10., 2.)";
      setMinGMean_ = false;
      setMaxGMean_ = true;
      fitParamCoeffGMean_.resize(2);
      momParamGMean_ = bookTF1("GMean", momParametrizationFormulaGMean_);
    } else if ( decayMode_ == kElectron_Muon && isSelected_ ) {
      momParametrizationFormulaGMean_ = 
	"[0] + [1]*TMath::Power(TMath::Max(1., x - [2]), [3]) + [4]*0.5*(1. + TMath::Erf([5]*(x - [6])))";
      setMinGMean_ = false;
      setMaxGMean_ = true;
      fitParamCoeffGMean_.resize(7);
      momParamGMean_ = bookTF1("GMean", momParametrizationFormulaGMean_);
      momParamGMean_->SetParameter(0,   3.);
      momParamGMean_->SetParameter(1,   1.);
      momParamGMean_->SetParameter(2,  10.);
      momParamGMean_->SetParameter(3,  -0.5);
      momParamGMean_->SetParameter(4,   0.5);
      momParamGMean_->SetParameter(5,   0.025);
      momParamGMean_->SetParameter(6,  75.);
    } else if ( (decayMode_ == kOneProng0Pi0   ||
		 decayMode_ == kOneProngGt0Pi0 ||
		 decayMode_ == kThreeProng0Pi0) && isSelected_ ) {
      momParametrizationFormulaGMean_ = "[0] + [1]/TMath::Power(x - 15., 2.)";
      setMinGMean_ = false;
      setMaxGMean_ = true;
      fitParamCoeffGMean_.resize(2);
      momParamGMean_ = bookTF1("GMean", momParametrizationFormulaGMean_);
    } else {
      momParametrizationFormulaGMean_ = "[0]";
      setMinGMean_ = false;
      setMaxGMean_ = false;
      fitParamCoeffGMean_.resize(1);
      momParamGMean_ = bookTF1("GMean", momParametrizationFormulaGMean_);
    }
    if ( decayMode_ == kElectron_Muon && isAll_ ) {
      momParametrizationFormulaGSigma_ = "[0]";
      setMinGSigma_ = false;
      setMaxGSigma_ = false;
      fitParamCoeffGSigma_.resize(1);
      momParamGSigma_ = bookTF1("GSigma", momParametrizationFormulaGSigma_);
    } else if ( decayMode_ == kElectron_Muon && isSelected_ ) {
      momParametrizationFormulaGSigma_ = 
	"[0] + [1]*TMath::Power(TMath::Max(1., x - [2]), [3]) + [4]*0.5*(1. + TMath::Erf([5]*(x - [6])))";
      setMinGSigma_ = false;
      setMaxGSigma_ = true;
      fitParamCoeffGSigma_.resize(7);
      momParamGSigma_ = bookTF1("GSigma", momParametrizationFormulaGSigma_);
      momParamGSigma_->SetParameter( 0,      1.    );
      momParamGSigma_->SetParameter( 1,      2.    );
      momParamGSigma_->SetParameter( 2,     10.    );
      momParamGSigma_->SetParameter( 3,     -1.    );
      momParamGSigma_->SetParameter( 4,      0.1   );
      momParamGSigma_->SetParameter( 5,      0.025 );
      momParamGSigma_->SetParameter( 6,     75.    );
    } else if ( decayMode_ == kOneProng0Pi0 && isSelected_ ) {
      momParametrizationFormulaGSigma_ = "[0] + [1]/TMath::Power(x - 15, 2.)";
      setMinGSigma_ = false;
      setMaxGSigma_ = true;
      fitParamCoeffGSigma_.resize(2);
      momParamGSigma_ = bookTF1("GSigma", momParametrizationFormulaGSigma_);
    } else if ( decayMode_ == kThreeProng0Pi0 && isSelected_ ) {
      momParametrizationFormulaGSigma_ = "[0] + [1]*x";
      setMinGSigma_ = false;
      setMaxGSigma_ = false;
      fitParamCoeffGSigma_.resize(2);
      momParamGSigma_ = bookTF1("GSigma", momParametrizationFormulaGSigma_);
    } else {
      momParametrizationFormulaGSigma_ = "([0] + [1]*x)*(1.0 + [2]*TMath::TanH([3] + [4]*x))";
      setMinGSigma_ = false;
      setMaxGSigma_ = false;
      fitParamCoeffGSigma_.resize(5);
      momParamGSigma_ = bookTF1("GSigma", momParametrizationFormulaGSigma_);
      momParamGSigma_->SetParameter( 0,      1.   );
      momParamGSigma_->SetParameter( 1,      0.   );
      momParamGSigma_->SetParameter( 2,     +0.2  );
      momParamGSigma_->SetParameter( 3, -60.*0.05 );
      momParamGSigma_->SetParameter( 4,      0.05 );
    }
    momParametrizationFormulaAlpha_ = "[0] + [1]*x";
    setMinAlpha_ = false;
    setMaxAlpha_ = false;
    fitParamCoeffAlpha_.resize(2);
    momParamAlpha_ = bookTF1("Alpha", momParametrizationFormulaAlpha_);
    if ( decayMode_ == kOneProng0Pi0 && isSelected_ ) {
      momParametrizationFormulaSlope_ = "[0]";
      setMinSlope_ = false;
      setMaxSlope_ = false;
      fitParamCoeffSlope_.resize(1);
      momParamSlope_ = bookTF1("Slope", momParametrizationFormulaSlope_);
    } else if ( decayMode_ == kOneProngGt0Pi0 && isSelected_ ) {
      momParametrizationFormulaSlope_ = "[0] + [1]/TMath::Power(x - 15., 2.)";
      setMinSlope_ = true;
      setMaxSlope_ = false;
      fitParamCoeffSlope_.resize(2);
      momParamSlope_ = bookTF1("Slope", momParametrizationFormulaSlope_);
    } else {
      momParametrizationFormulaSlope_ = "[0] + [1]*x";
      setMinSlope_ = false;
      setMaxSlope_ = false;
      fitParamCoeffSlope_.resize(2);
      momParamSlope_ = bookTF1("Slope", momParametrizationFormulaSlope_);
    }
    if ( decayMode_ == kOneProng0Pi0 && isSelected_ ) {
      momParametrizationFormulaOffset_ = "[0] + [1]/TMath::Power(x - 15., 2.)";
      setMinOffset_ = false;
      setMaxOffset_ = true;
      fitParamCoeffOffset_.resize(2);
      momParamOffset_ = bookTF1("Offset", momParametrizationFormulaOffset_);
    } else {
      momParametrizationFormulaOffset_ = "[0] + [1]*x";
      setMinOffset_ = false;
      setMaxOffset_ = false;
      fitParamCoeffOffset_.resize(2);
      momParamOffset_ = bookTF1("Offset", momParametrizationFormulaOffset_);
    }
    momParametrizationFormulaC_       = 
      "0.5*(1.0 + TMath::Cos([0])) + 0.125*(1.0 - TMath::Cos([0]))*(1.0 + TMath::Cos([1]))*(1.0 + TMath::TanH([2] + [3]*x))";
    setMinC_ = false;
    setMaxC_ = false;
    fitParamCoeffC_.resize(4);
    momParamC_ = bookTF1("C", momParametrizationFormulaC_);
    momParamC_->SetParameter( 0,      0.5*TMath::Pi() );
    momParamC_->SetParameter( 1,      0.5*TMath::Pi() );
    momParamC_->SetParameter( 2, -60.*0.05            );
    momParamC_->SetParameter( 3,      0.05            );
//--- parametrize pt/energy dependence of MP1 fitParameter by (orthogonal) Chebyshev polynomials
//   ( cf. http://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html )
    momParametrizationFormulaMP2_     =
      "TMath::Abs((1./(x*x*x*x))*([0] + [1]*x + [2]*(2.0*x*x - 1.0) + [3]*(4.0*x*x*x - 3.0*x) + [4]*(8.0*x*x*x*x - 8.0*x*x + 1.0)))";
    setMinMP2_ = false;
    setMaxMP2_ = false;
    fitParamCoeffMP2_.resize(5);
    momParamMP2_ = bookTF1("MP2", momParametrizationFormulaMP2_);
    if ( decayMode_ == kElectron_Muon && isAll_ ) {
      momParametrizationFormulaWidth2_ = "([0] - [1])*0.5*(1. + TMath::Erf([2]*(x - [3]))) + [1]";
      setMinWidth2_ = false;
      setMaxWidth2_ = false;
      fitParamCoeffWidth2_.resize(4);
      momParamWidth2_ = bookTF1("Width2", momParametrizationFormulaWidth2_);
      momParamWidth2_->SetParameter( 0,    1.e-1 );
      momParamWidth2_->SetParameter( 1,    1.e-2 );
      momParamWidth2_->SetParameter( 2,    1.e-2 );
      momParamWidth2_->SetParameter( 3,   75. );
    } else if ( decayMode_ == kElectron_Muon && isSelected_ ) {
      momParametrizationFormulaWidth2_ = 
	"[0] + [1]/TMath::Power(x - [2], [3]) + [4]*0.5*(1. + TMath::Erf([5]*(x - [6])))";
      setMinWidth2_ = false;
      setMaxWidth2_ = true;
      fitParamCoeffWidth2_.resize(7);
      momParamWidth2_ = bookTF1("Width2", momParametrizationFormulaWidth2_);
    } else if ( decayMode_ == kOneProng0Pi0 && isSelected_ ) {      
      momParametrizationFormulaWidth2_ = "[0] + [1]/TMath::Power(x - 15., 2.)";
      setMinWidth2_ = false;
      setMaxWidth2_ = true;
      fitParamCoeffWidth2_.resize(2);
      momParamWidth2_ = bookTF1("Width2", momParametrizationFormulaWidth2_);
    } else if ( decayMode_ == kOneProngGt0Pi0 && isSelected_ ) {      
      momParametrizationFormulaWidth2_ = "[0]";
      setMinWidth2_ = false;
      setMaxWidth2_ = false;
      fitParamCoeffWidth2_.resize(1);
      momParamWidth2_ = bookTF1("Width2", momParametrizationFormulaWidth2_);
    } else if ( decayMode_ == kThreeProng0Pi0 && isSelected_ ) {
      momParametrizationFormulaWidth2_ = "[0] + [1]/TMath::Power(x - 15., 2.)";
      setMinWidth2_ = false;
      setMaxWidth2_ = true;
      fitParamCoeffWidth2_.resize(2);
      momParamWidth2_ = bookTF1("Width2", momParametrizationFormulaWidth2_);
    } else {
      momParametrizationFormulaWidth2_ = "TMath::Max(1.e-4, [0] + [1]*x)";
      setMinWidth2_ = false;
      setMaxWidth2_ = false;
      fitParamCoeffWidth2_.resize(2);
      momParamWidth2_ = bookTF1("Width2", momParametrizationFormulaWidth2_);
    }
    if ( decayMode_ == kElectron_Muon && isAll_ ) {
      momParametrizationFormulaMP3_ = "(1. + [4]*x)*(([0] - [1])*0.5*(1. + TMath::Erf([2]*(x - [3]))) + [1])";
      setMinMP3_ = false;
      setMaxMP3_ = false;
      fitParamCoeffMP3_.resize(5);
      momParamMP3_ = bookTF1("MP3", momParametrizationFormulaMP3_);
      momParamMP3_->SetParameter( 0,       1.e-3 );
      momParamMP3_->SetParameter( 1,       1.e-4 );
      momParamMP3_->SetParameter( 2,      -1.e-2 );
      momParamMP3_->SetParameter( 3,      -2.e+2 );
      momParamMP3_->SetParameter( 4,       2.    );
    } else if ( decayMode_ == kElectron_Muon && isSelected_ ) {
      momParametrizationFormulaMP3_ = "([0] + [1]*x)*(1.0 + [2]*TMath::TanH([3] + [4]*x))";
      setMinMP3_ = false;
      setMaxMP3_ = false;
      fitParamCoeffMP3_.resize(5);
      momParamMP3_ = bookTF1("MP3", momParametrizationFormulaMP3_);
      momParamMP3_->SetParameter( 0,       7.5  );
      momParamMP3_->SetParameter( 1,       0.   );
      momParamMP3_->SetParameter( 2,      -0.3  );
      momParamMP3_->SetParameter( 3, -125.*0.10 );
      momParamMP3_->SetParameter( 4,       0.10 );
    } else {
      //momParametrizationFormulaMP3_ = "[0]*TMath::Min(1., [1]*(x - [2]))";
      //fitParamCoeffMP3_.resize(3);
      momParametrizationFormulaMP3_ = "([0] - [1])*0.5*(1. + TMath::Erf([2]*(x - [3]))) + [1]";
      setMinMP3_ = false;
      setMaxMP3_ = false;
      fitParamCoeffMP3_.resize(4);
      momParamMP3_ = bookTF1("MP3", momParametrizationFormulaMP3_);
      //momParamMP3_->SetParameter( 0,  9.5  );
      //momParamMP3_->SetParameter( 1,  0.01 );
      //momParamMP3_->SetParameter( 2, 10.   );
      momParamMP3_->SetParameter( 0, 10.    );
      momParamMP3_->SetParameter( 1,  0.    );
      momParamMP3_->SetParameter( 2,  0.025 );
      momParamMP3_->SetParameter( 3, 50.    );
    }
    if ( decayMode_ == kElectron_Muon && isSelected_ ) {
      momParametrizationFormulaWidth3_ = "([0] + [1]*x)*(1.0 + [2]*TMath::TanH([3] + [4]*x))";
      setMinWidth3_ = false;
      setMaxWidth3_ = false;
      fitParamCoeffWidth3_.resize(5);
      momParamWidth3_ = bookTF1("Width3", momParametrizationFormulaWidth3_);
      momParamWidth3_->SetParameter( 0,       0.1  );
      momParamWidth3_->SetParameter( 1,       0.   );
      momParamWidth3_->SetParameter( 2,      +0.5  );
      momParamWidth3_->SetParameter( 3, -250.*0.10 );
      momParamWidth3_->SetParameter( 4,       0.10 );
    } else if ( decayMode_ == kOneProngGt0Pi0 && isSelected_ ) {
      momParametrizationFormulaWidth3_ = "[0]";
      setMinWidth3_ = false;
      setMaxWidth3_ = false;
      fitParamCoeffWidth3_.resize(1);
      momParamWidth3_ = bookTF1("Width3", momParametrizationFormulaWidth3_);
    } else {
      momParametrizationFormulaWidth3_ = "TMath::Max(5.e-3, [0] + [1]*x)";
      setMinWidth3_ = false;
      setMaxWidth3_ = false;
      fitParamCoeffWidth3_.resize(2);
      momParamWidth3_ = bookTF1("Width3", momParametrizationFormulaWidth3_);
    }
    if ( decayMode_ == kOneProng0Pi0 && isSelected_ ) {
      momParametrizationFormulaX0_ = "[0]*0.5*(1. + TMath::Erf([1]*(x - [2])))";
      setMinX0_ = false;
      setMaxX0_ = false;
      fitParamCoeffX0_.resize(3);
      momParamX0_ = bookTF1("X0", momParametrizationFormulaX0_);
      momParamX0_->SetParameter( 0,       3.  );
      momParamX0_->SetParameter( 1,       1.  );
      momParamX0_->SetParameter( 2,      30.  );
    } else {
      momParametrizationFormulaX0_ = "[0] + [1]*x";
      setMinX0_ = false;
      setMaxX0_ = false;
      fitParamCoeffX0_.resize(2);
      momParamX0_ = bookTF1("X0", momParametrizationFormulaX0_);
    }
    if ( decayMode_ == kElectron_Muon && isSelected_ ) {
//--- parametrize pt/energy dependence of X0 fitParameter by (orthogonal) Chebyshev polynomials
//   ( cf. http://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html )
      momParametrizationFormulaDeltaX1_ = 
	"TMath::Abs((1./(x*x*x*x))*([0] + [1]*x + [2]*(2.0*x*x - 1.0) + [3]*(4.0*x*x*x - 3.0*x) + [4]*(8.0*x*x*x*x - 8.0*x*x + 1.0)))";
      setMinDeltaX1_ = false;
      setMaxDeltaX1_ = false;
      fitParamCoeffDeltaX1_.resize(5);
      momParamDeltaX1_ = bookTF1("DeltaX1", momParametrizationFormulaDeltaX1_);
      momParamDeltaX1_->SetParameter( 0,       1.   );
      momParamDeltaX1_->SetParameter( 1,       0.   );
      momParamDeltaX1_->SetParameter( 2,       2.5  );
      momParamDeltaX1_->SetParameter( 3, -125.*0.10 );
      momParamDeltaX1_->SetParameter( 4,       0.10 );
    } else if ( decayMode_ == kOneProng0Pi0 && isSelected_ ) {
      momParametrizationFormulaDeltaX1_ = "[0] + [1]/TMath::Power(x - 15., 2.)";
      setMinDeltaX1_ = false;
      setMaxDeltaX1_ = true;
      fitParamCoeffDeltaX1_.resize(2);
      momParamDeltaX1_ = bookTF1("DeltaX1", momParametrizationFormulaDeltaX1_);
    } else if ( (decayMode_ == kOneProngGt0Pi0 ||
		 decayMode_ == kThreeProng0Pi0) && isSelected_ ) {
      momParametrizationFormulaDeltaX1_ = "[0] + [1]/TMath::Power(x - 15., 2.)";
      setMinDeltaX1_ = false;
      setMaxDeltaX1_ = true;
      fitParamCoeffDeltaX1_.resize(2);
      momParamDeltaX1_ = bookTF1("DeltaX1", momParametrizationFormulaDeltaX1_);
    } else {
      momParametrizationFormulaDeltaX1_ = "[0] + [1]*x";
      setMinDeltaX1_ = false;
      setMaxDeltaX1_ = false;
      fitParamCoeffDeltaX1_.resize(2);
      momParamDeltaX1_ = bookTF1("DeltaX1", momParametrizationFormulaDeltaX1_);
    }
    if ( decayMode_ == kElectron_Muon && isSelected_ ) {
      momParametrizationFormulaDeltaX2_ = "([0] + [1]*x)*(1.0 + [2]*TMath::TanH([3] + [4]*x))";
      setMinDeltaX2_ = false;
      setMaxDeltaX2_ = false;
      fitParamCoeffDeltaX2_.resize(5);
      momParamDeltaX2_ = bookTF1("DeltaX2", momParametrizationFormulaDeltaX2_);
      momParamDeltaX2_->SetParameter( 0,       5.   );
      momParamDeltaX2_->SetParameter( 1,       0.   );
      momParamDeltaX2_->SetParameter( 2,      -0.5  );
      momParamDeltaX2_->SetParameter( 3, -125.*0.10 );
      momParamDeltaX2_->SetParameter( 4,       0.10 );
    } else {
      momParametrizationFormulaDeltaX2_ = "([0] - [1])*0.5*(1. + TMath::Erf([2]*(x - [3]))) + [1]";
      if ( decayMode_ == kOneProngGt0Pi0 && isSelected_ ) setMinDeltaX2_ = true;
      else setMinDeltaX2_ = false;
      setMaxDeltaX2_ = false;
      fitParamCoeffDeltaX2_.resize(4);
      momParamDeltaX2_ = bookTF1("DeltaX2", momParametrizationFormulaDeltaX2_);
      momParamDeltaX2_->SetParameter( 0, 10.    );
      momParamDeltaX2_->SetParameter( 1,  0.    );
      momParamDeltaX2_->SetParameter( 2,  0.025 );
      momParamDeltaX2_->SetParameter( 3, 50.    );
    }
  }

  ~fitManager()
  {
    clearCollection(prefitParamMP1_);
    clearCollection(prefitParamWidth1_);
    clearCollection(prefitParamGMean_);
    clearCollection(prefitParamGSigma_);
    clearCollection(prefitParamAlpha_);
    clearCollection(prefitParamSlope_);
    clearCollection(prefitParamOffset_);
    clearCollection(prefitParamC_);
    clearCollection(prefitParamMP2_);
    clearCollection(prefitParamWidth2_);
    clearCollection(prefitParamMP3_);
    clearCollection(prefitParamWidth3_);
    clearCollection(prefitParamX0_);
    clearCollection(prefitParamDeltaX1_);    
    clearCollection(prefitParamDeltaX2_);   
    clearCollection(prefitModel_);

    delete prefitResultMP1_;
    delete prefitResultWidth1_;
    delete prefitResultGMean_;
    delete prefitResultGSigma_;
    delete prefitResultAlpha_;
    delete prefitResultSlope_;
    delete prefitResultOffset_;
    delete prefitResultC_;
    delete prefitResultMP2_;
    delete prefitResultWidth2_;
    delete prefitResultMP3_;
    delete prefitResultWidth3_;
    delete prefitResultX0_;
    delete prefitResultDeltaX1_;
    delete prefitResultDeltaX2_;

    delete momParamMP1_;
    delete momParamWidth1_;
    delete momParamGMean_;
    delete momParamGSigma_;
    delete momParamAlpha_;
    delete momParamSlope_;
    delete momParamOffset_;
    delete momParamC_;
    delete momParamMP2_;
    delete momParamWidth2_;
    delete momParamMP3_;
    delete momParamWidth3_;
    delete momParamX0_;
    delete momParamDeltaX1_;
    delete momParamDeltaX2_;

    clearCollection(fitParamCoeffMP1_);
    delete fitParamMP1_;
    clearCollection(fitParamCoeffWidth1_);
    delete fitParamWidth1_;
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
    clearCollection(fitParamCoeffMP2_);
    delete fitParamMP2_;
    clearCollection(fitParamCoeffWidth2_);
    delete fitParamWidth2_;
    clearCollection(fitParamCoeffMP3_);
    delete fitParamMP3_;
    clearCollection(fitParamCoeffWidth3_);
    delete fitParamWidth3_;
    clearCollection(fitParamCoeffX0_);
    delete fitParamX0_;
    clearCollection(fitParamCoeffDeltaX1_);
    delete fitParamDeltaX1_;
    clearCollection(fitParamCoeffDeltaX2_);
    delete fitParamDeltaX2_;
    delete fitModel_;
  }

  TF1* bookTF1(const TString& fitParameter, const TString& formula)
  {
    std::cout << "<fitManager::bookTF1>:" << std::endl;
    std::cout << " " << fitParameter.Data() << ": formula = " << formula.Data() << std::endl;

    TString decayMode_string = getDecayMode_string(decayMode_);
    TString tf1Name = Form("%s_%s_%s_%s", fitParameter.Data(), sep_->GetName(), decayMode_string.Data(), label_.Data());
    double momMin = momBinning_[0];
    //std::cout << "momMin = " << momMin << std::endl;
    double momMax = momBinning_[momBinning_.GetSize() - 1];
    //std::cout << "momMax = " << momMax << std::endl;
    TF1* tf1 = new TF1(tf1Name.Data(), formula, momMin, momMax);
    return tf1;
  }

  void fitTF1(TF1* tf1, bool limitMinValue, bool limitMaxValue, TGraph* graph, 
	      std::vector<double>& fitParamValues, std::vector<double>& fitParamErrors)
  {
    std::cout << "<fitManager::fitTF1>:" << std::endl;
    std::cout << " " << tf1->GetName() << ": formula = " << tf1->GetTitle() << std::endl;

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

    // CV: if TF1 formula contains reserved keywords 'MIN'/'MAX'
    //     replace 'MIN'/'MAX' keywords by minimum/maximum y-value of TGraph object
    //    (this feature prevents extrapolation of TF1 beyond fitted range,
    //     which may be important in particular when fitting functions like ~1./x^N, N > 0.)
    double minValue = 0.;
    double maxValue = 0.;

    int numPoints = graph->GetN();
    for ( int iPoint = 0; iPoint < numPoints; ++iPoint ) {
      double x, y;
      graph->GetPoint(iPoint, x, y);
      
      if ( iPoint == 0 ) {
	minValue = y;
	maxValue = y;
      } else {
	if ( y < minValue ) minValue = y;
	if ( y > maxValue ) maxValue = y;
      }
    }

    std::cout << "minValue = " << minValue << std::endl;
    std::cout << "maxValue = " << maxValue << std::endl;

    if ( limitMinValue ) {
      TString formula_old = tf1->GetTitle();
      TString formula_new = Form("TMath::Max(%e, %s)", minValue, formula_old.Data());
      tf1->SetTitle(formula_new.Data());
      std::cout << "--> replacing formula = " << tf1->GetTitle() << std::endl;
      tf1->Compile();
    }

    if ( limitMaxValue ) {
      TString formula_old = tf1->GetTitle();
      TString formula_new = Form("TMath::Min(%e, %s)", maxValue, formula_old.Data());
      tf1->SetTitle(formula_new.Data());
      std::cout << "--> replacing formula = " << tf1->GetTitle() << std::endl;
      tf1->Compile();
    }

    for ( int iParameter = 0; iParameter < numFitParameter; ++iParameter ) {
      tf1->SetParameter(iParameter, fitParamValues[iParameter]);
      tf1->SetParError(iParameter,  fitParamErrors[iParameter]);
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

  RooRealVar* buildPrefitParameter(const TString& name, const TString& momBinName, double startValue, 
				   double min, double max, bool isConstant)
  {
    std::cout << "<fitManager::buildPrefitParameter>:" << std::endl;

    TString decayMode_string = getDecayMode_string(decayMode_);
    TString fitParameterName = 
      Form("%s_%s_%s_%s_%s", name.Data(), momBinName.Data(), sep_->GetName(), decayMode_string.Data(), label_.Data());
    RooRealVar* fitParameter = new RooRealVar(fitParameterName.Data(), fitParameterName.Data(), startValue, min, max);
    if ( isConstant ) fitParameter->setConstant();
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
    printPrefitParameter("mp1", prefitParamMP1_[momBin]);
    printPrefitParameter("width1", prefitParamWidth1_[momBin]);
    printPrefitParameter("gmean", prefitParamGMean_[momBin]);
    printPrefitParameter("gsigma", prefitParamGSigma_[momBin]);
    printPrefitParameter("alpha", prefitParamAlpha_[momBin]);
    printPrefitParameter("slope", prefitParamSlope_[momBin]);
    printPrefitParameter("offset", prefitParamOffset_[momBin]);
    printPrefitParameter("C", prefitParamC_[momBin]);
    printPrefitParameter("mp2", prefitParamMP2_[momBin]);
    printPrefitParameter("width2", prefitParamWidth2_[momBin]);
    printPrefitParameter("mp3", prefitParamMP3_[momBin]);
    printPrefitParameter("width3", prefitParamWidth3_[momBin]);
    printPrefitParameter("x0", prefitParamX0_[momBin]);
    printPrefitParameter("dx1", prefitParamDeltaX1_[momBin]);
    printPrefitParameter("dx2", prefitParamDeltaX2_[momBin]);
  }

  void runPrefit(const TString& inputFileName, const TString& inputDirName, const TString& label, const TString& outputFileName)
  {
    std::cout << "<fitManager::runPrefit>:" << std::endl;
    std::cout << " inputFileName = " << inputFileName.Data() << std::endl;
    std::cout << " inputDirName = " << inputDirName.Data() << std::endl;
    
    printTimeStamp("<fitManager::runPrefit>");

    TString decayMode_string = getDecayMode_string(decayMode_);

    RooRealVar* sepTimesMom_prefit = new RooRealVar("sepTimesMom_prefit", "sepTimesMom_prefit", 0., 25.);
    sepTimesMom_prefit->setBins(250);

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
      std::cout << " histogramMean = " << histogramMean << std::endl;
      double histogramMax_x = histogram->GetBinCenter(histogram->GetMaximumBin());
      std::cout << " histogramMax_x = " << histogramMax_x << std::endl;
      double histogramMax_y = histogram->GetBinContent(histogram->GetMaximumBin());
      std::cout << " histogramMax_y = " << histogramMax_y << std::endl;
      double histogramMax_binWidth = histogram->GetBinWidth(histogram->GetMaximumBin());
      double errEdgePosRight = 0.;
      double errEdgePosLeft  = 0.;
      double fallingEdge_position = compFallingEdgePos(histogram, decayMode_, errEdgePosRight, errEdgePosLeft);
      std::cout << " fallingEdge_position = " << fallingEdge_position << std::endl;    
      double prePeak_x_left  = 0.;
      double prePeak_x_right = 0.;
      compPrePeakPos(histogram, decayMode_, prePeak_x_left, prePeak_x_right);
      std::cout << " prePeak_x = " << prePeak_x_left << ".." << prePeak_x_right << std::endl;  
      int prePeak_bin_right = histogram->FindBin(prePeak_x_right);
      double prePeak_y_right = histogram->GetBinContent(prePeak_bin_right);
      std::cout << " prePeak_y(right) = " << prePeak_y_right << std::endl; 
      bool prePreak_exists = (prePeak_x_right >= (prePeak_x_left + 0.1));
      double histogramRMS = histogram->GetRMS();
      //if ( histogramRMS  > 0.15 ) histogramRMS  = 0.15;
      std::cout << " histogramRMS = " << histogramRMS << std::endl;
      double histogramRMSltMax = compHistogramRMSltMax(histogram);
      std::cout << " histogramRMSltMax = " << histogramRMSltMax << std::endl;
      double histogramRMSgtMax = compHistogramRMSgtMax(histogram);
      std::cout << " histogramRMSgtMax = " << histogramRMSgtMax << std::endl;
      double histogramIntegral = histogram->Integral();
      std::cout << " histogramIntegral = " << histogramIntegral << std::endl;

      const double epsilon = 1.e-3;

      //-------------------------------------------------------------------------
      // set parameters of Landau1
      double mp1Initial = -1.;
      double mp1Min     = mp1Initial - epsilon;
      double mp1Max     = mp1Initial + epsilon;
      double mp1_isConstant = true;
      if ( decayMode_ == kOneProng0Pi0 && prePreak_exists ) {
	mp1Initial = prePeak_x_left;
	mp1Min     = mp1Initial - 2.*histogram->GetBinWidth(histogram->GetMaximumBin());
        mp1Max     = TMath::Max(prePeak_x_left, prePeak_x_right - 1.*histogram->GetBinWidth(histogram->GetMaximumBin()));
	mp1_isConstant = false;
      }
      prefitParamMP1_[iMomBin] = buildPrefitParameter("mp1", momBinName, mp1Initial, mp1Min, mp1Max, mp1_isConstant);
      double width1Initial = 0.25;
      double width1Min     = width1Initial - epsilon;
      double width1Max     = width1Initial + epsilon;
      bool width1_isConstant = true;
      if ( decayMode_ == kOneProng0Pi0 && prePreak_exists ) {
	width1Initial = 0.25;
	width1Min     = 0.15;
	width1Max     = 0.50;
	width1_isConstant = false;
      }
      prefitParamWidth1_[iMomBin] = buildPrefitParameter("width1", momBinName, width1Initial, width1Min, width1Max, width1_isConstant);
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // set parameters of Gaussian
      double gmeanInitial = histogramMax_x;
      double gmeanMin     = 0.5*histogramMax_x;
      double gmeanMax     = 1.5*histogramMax_x;
      bool gmean_isConstant = false;
      if ( decayMode_ == kOneProng0Pi0 ) {
	gmeanInitial = histogramMax_x;
	gmeanMin     = histogramMax_x - 2.*histogramMax_binWidth;
	gmeanMax     = histogramMax_x + 2.*histogramMax_binWidth;
      }
      prefitParamGMean_[iMomBin] = buildPrefitParameter("gmean", momBinName, gmeanInitial, gmeanMin, gmeanMax, gmean_isConstant);
      double gsigmaInitial = histogramRMS; 
      double gsigmaMin     = 0.25*TMath::Min(histogramRMSltMax, histogramRMS);
      double gsigmaMax     = 2.0*TMath::Max(histogramRMSltMax, histogramRMS);
      bool gsigma_isConstant = false;
      if ( decayMode_ == kOneProng0Pi0 ) {
	gsigmaInitial = histogramMax_x - prePeak_x_right;
	gsigmaMin     = 2.5e-2;
	gsigmaMax     = 2.*gsigmaInitial;
      }
      prefitParamGSigma_[iMomBin] = buildPrefitParameter("gsigma", momBinName, gsigmaInitial, gsigmaMin, gsigmaMax, gsigma_isConstant);
      double alphaInitial =  0.;
      double alphaMin     = -1.e-9;
      double alphaMax     = +1.e-9;
      bool alpha_isConstant = true;
      //if ( decayMode_ == kElectron_Muon ) {
      //  alphaInitial = -5.;
      //  alphaMin     = -1.e+1;
      //  alphaMax     = +1.e+1;
      //  alpha_isConstant = false;	
      if ( decayMode_ == kOneProngGt0Pi0 && isAll_ ) {
	alphaInitial =  0.;
	alphaMin     = -1.e-1;
	alphaMax     = +1.e+1;
	alpha_isConstant = false;
      } else if ( decayMode_ == kThreeProng0Pi0 ) {
	alphaInitial = -5.;
	alphaMin     = -1.e+1;
	alphaMax     = +1.e+1;
	alpha_isConstant = false;
      }
      prefitParamAlpha_[iMomBin] = buildPrefitParameter("alpha", momBinName, alphaInitial, alphaMin, alphaMax, alpha_isConstant);
      if ( alpha_isConstant ) prefitParamAlpha_[iMomBin]->setConstant();
      double slopeInitial = 1./histogramMax_x;
      double slopeMin     = 1.e-2*slopeInitial;
      double slopeMax     = 1.e+2*slopeInitial;
      bool slope_isConstant = false;
      if ( decayMode_ == kElectron_Muon ) {
	slopeInitial = 1./histogramMax_x;
        slopeMin     = 0.1*slopeInitial;
	slopeMax     = 2.*slopeInitial;
      } else if ( decayMode_ == kOneProng0Pi0 ) {
	if ( prePreak_exists ) {
	  //slopeInitial = (1. - prePeak_y_right/histogramMax_y)/(histogramMax_x - prePeak_x_right);
	  //slopeMin     = 0.5*slopeInitial;
	  //slopeMax     = 2.*slopeInitial;
	  slopeInitial = 1./histogramRMSltMax;
	  slopeMin     = 1.e-2*slopeInitial;
	  slopeMax     = 1.e+2*slopeInitial;
	} else {
	  slopeInitial = 1./histogramRMSltMax;
	  slopeMin     = 1.e-2*slopeInitial;
	  slopeMax     = 1.e+2*slopeInitial;
	}
      } else if ( decayMode_ == kOneProngGt0Pi0 && isSelected_ ) {
	slopeInitial = 1./histogramMax_x;
        slopeMin     = 0.1*slopeInitial;
	slopeMax     = 2.*slopeInitial;
      } else if ( decayMode_ == kOneProngGt0Pi0 ||
		  decayMode_ == kThreeProng0Pi0 ) {
	slopeInitial =  0.;
	slopeMin     = -1.e-9;
	slopeMax     = +1.e-9;
	slope_isConstant = true;
      } 
      prefitParamSlope_[iMomBin] = buildPrefitParameter("slope", momBinName, slopeInitial, slopeMin, slopeMax, slope_isConstant);
      double offsetInitial =  0.;
      double offsetMin     = -1.e-9;
      double offsetMax     = +1.e-9;
      bool offset_isConstant = true;
      if ( slopeInitial > 0. ) {
	if ( prePreak_exists ) {
	  offsetInitial = prePeak_x_right;
	  offsetMin     = prePeak_x_right - 1.*histogramMax_binWidth;
	  offsetMax     = prePeak_x_right + 1.*histogramMax_binWidth;
	} else {
	  //offsetInitial = histogramMax_x - histogramMax_y/slopeInitial;
	  offsetInitial = histogramMax_x - 3.*histogramRMSltMax;
	  offsetMin     = 0.;
	  offsetMax     = histogramMax_x;
	}
	offset_isConstant = (gmean_isConstant || slope_isConstant);
      }
      prefitParamOffset_[iMomBin] = buildPrefitParameter("offset", momBinName, offsetInitial, offsetMin, offsetMax, offset_isConstant);
      double Cinitial = 0.25;
      double Cmin     = 0.;
      double Cmax     = 1.;
      bool C_isConstant = false;
      if ( decayMode_ == kElectron_Muon ) {
	if ( isAll_ ) {
	  Cinitial = 0.5;
	  Cmin     = 0.25;
	  Cmax     = 1.;
	} else if ( isSelected_ ) {
	  Cinitial = 1.;
	  Cmin     = Cinitial - epsilon;
	  Cmax     = Cinitial + epsilon;
	  C_isConstant = true;
	  //Cinitial = 0.5;
	  //Cmin     = 0.25;
	  //Cmax     = 1.;
	}
      } else if ( decayMode_ == kOneProng0Pi0 ) {
	//if ( (histogramMax_x - prePeak_x_right) >= 2.*histogramMax_binWidth ) {
	//  Cinitial          = 1.;
	//  Cmin              = Cinitial - epsilon;
	//  Cmax              = Cinitial + epsilon;
	//  slope_isConstant  = true;
	//  offset_isConstant = true;
	//} else {
	  Cinitial          = 0.;
	  Cmin              = Cinitial - epsilon;
	  Cmax              = Cinitial + epsilon;
	  gmean_isConstant  = true;
	  gsigma_isConstant = true;
	//}
	C_isConstant = true;
      } else if ( decayMode_ == kOneProngGt0Pi0 ||
		  decayMode_ == kThreeProng0Pi0 ) {
        Cinitial = 1.;
	Cmin     = 1. - epsilon;
	Cmax     = 1. + epsilon;
	C_isConstant = true;
      }
      prefitParamC_[iMomBin] = buildPrefitParameter("C", momBinName, Cinitial, Cmin, Cmax, C_isConstant);
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // set parameters of Landau2
      double mp2Initial = histogramMax_x;
      double mp2Min     = 0.8*histogramMax_x;
      double mp2Max     = fallingEdge_position;
      double mp2_isConstant = false;
      if ( decayMode_ == kOneProng0Pi0 && isAll_ ) {
	mp2Initial = histogramMax_x - 1.;
	mp2Min     = mp2Initial - 1.;
	mp2Max     = mp2Initial + 1.;
      } else if ( decayMode_ == kOneProng0Pi0 && isSelected_ ) {
	mp2Initial = histogramMax_x;
	mp2Min     = mp2Initial - 1.*histogramMax_binWidth;
	mp2Max     = mp2Initial + 1.*histogramMax_binWidth;
	//mp2_isConstant = true;
      } else if ( decayMode_ == kOneProngGt0Pi0 && isSelected_ && isDeltaR_ ) {
	mp2Initial = 2.3;
	mp2Min     = 2.25;
	mp2Max     = 2.35;
      } else if ( decayMode_ == kThreeProng0Pi0 ) {
	mp2Initial = 0.5*(mp2Min + mp2Max);
	mp2Min     = 0.6*histogramMax_x;
	mp2Max     = 0.9*histogramMax_x;
      }
      prefitParamMP2_[iMomBin] = buildPrefitParameter("mp2", momBinName, mp2Initial, mp2Min, mp2Max, mp2_isConstant);
      double width2Initial = 0.04*histogramRMSgtMax;
      double width2Min     = 5.e-3;
      double width2Max     = 2.0*histogramRMSgtMax;
      bool width2_isConstant = false;
      if ( decayMode_ == kElectron_Muon ) {
	//width2Initial = 2.e-3;
	width2Min = 1.e-4;
	width2_isConstant = true;
      } else if ( decayMode_ == kOneProng0Pi0 && isAll_ ) {
	width2Initial = 1.e-1;
	width2Min     = 1.e-4;	
      } else if ( decayMode_ == kOneProng0Pi0 && isSelected_ ) {
	//width2Initial = 1.0*histogramRMSgtMax;
	width2Initial = 0.5;
	width2Min     = 1.e-4;
      } else if ( decayMode_ == kThreeProng0Pi0 ) {
	width2Min = 5.e-3;
	width2Max = 2.0*histogramRMSgtMax;
      }
      prefitParamWidth2_[iMomBin] = buildPrefitParameter("width2", momBinName, width2Initial, width2Min, width2Max, width2_isConstant);
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // set parameters of Landau3
      double mp3Initial = histogramMax_x + 2.;
      double mp3Min     = fallingEdge_position + 2.;
      double mp3Max     = 25.;
      double mp3_isConstant = false;
      if ( decayMode_ == kOneProng0Pi0 ) {
	mp3Initial = 3.;
        mp3Min     = 2.5;
        mp3Max     = 3.5;
      }
      prefitParamMP3_[iMomBin] = buildPrefitParameter("mp3", momBinName, mp3Initial, mp3Min, mp3Max, mp3_isConstant);
      double width3Initial = 5.e-1;
      double width3Min     = 5.e-3;
      double width3Max     = 2.e-0;
      bool width3_isConstant = false;
      prefitParamWidth3_[iMomBin] = buildPrefitParameter("width3", momBinName, width3Initial, width3Min, width3Max, width3_isConstant);
      //-------------------------------------------------------------------------

      double x0Initial =  0.;
      double x0Min     = -1.e-9;
      double x0Max     = +1.e-9;
      bool x0_isConstant = true;
      if ( decayMode_ == kOneProng0Pi0 && prePreak_exists ) {
	x0Initial = prePeak_x_right;
        x0Min     = x0Initial - 2.*histogramMax_binWidth;
        x0Max     = x0Initial + 2.*histogramMax_binWidth;
        x0_isConstant = false;
      }
      prefitParamX0_[iMomBin] = buildPrefitParameter("x0", momBinName, x0Initial, x0Min, x0Max, x0_isConstant);
      double dx1Initial = fallingEdge_position - x0Initial;
      double dx1Min     = dx1Initial - (0.1*fallingEdge_position + errEdgePosLeft);
      double dx1Max     = dx1Initial + (0.1*fallingEdge_position + errEdgePosRight);
      bool dx1_isConstant = false;
      if ( decayMode_ == kOneProng0Pi0 ) {
	if ( isAll_ ) {
	  dx1Initial = 0.5;
	  dx1Min     = dx1Initial - 0.25;
	  dx1Max     = dx1Initial + 0.35;
	} else if ( isSelected_ ) {
	  dx1Initial = histogramMax_x - x0Initial + 1.*histogramMax_binWidth;
	  dx1Min     = dx1Initial - 1.e-9;
	  dx1Max     = dx1Initial + 0.25;
	} else assert(0);
	dx1_isConstant = false;
      } else if ( decayMode_ == kThreeProng0Pi0 ) {
	dx1Min = dx1Initial - (0.2*fallingEdge_position + errEdgePosLeft);
	dx1Max = dx1Initial + (0.2*fallingEdge_position + errEdgePosRight);
      }
      prefitParamDeltaX1_[iMomBin] = buildPrefitParameter("dx1", momBinName, dx1Initial, dx1Min, dx1Max, dx1_isConstant);
      double dx2Initial = 25.;
      double dx2Min     = 25. - epsilon;
      double dx2Max     = 25. + epsilon;
      //bool dx2_isConstant = true;
      bool dx2_isConstant = false;
      if ( decayMode_ == kElectron_Muon && isAll_ ) {
	dx2Min         = 0.5;
      } else if ( decayMode_ == kElectron_Muon && isSelected_ ) {
	dx2Min         = 0.5;
	//dx2_isConstant = true;
      } else if ( decayMode_ == kOneProng0Pi0 ) {
	dx2_isConstant = true;
      } else if ( decayMode_ == kOneProngGt0Pi0 && isSelected_ && isDeltaR_ ) {
	dx2_isConstant = true;
      } 
      prefitParamDeltaX2_[iMomBin] = buildPrefitParameter("dx2", momBinName, dx2Initial, dx2Min, dx2Max, dx2_isConstant);

      prefitParamMP1_[iMomBin]->setConstant(mp1_isConstant);
      prefitParamWidth1_[iMomBin]->setConstant(width1_isConstant);
      prefitParamGMean_[iMomBin]->setConstant(gmean_isConstant);
      prefitParamGSigma_[iMomBin]->setConstant(gsigma_isConstant);
      prefitParamAlpha_[iMomBin]->setConstant(alpha_isConstant);
      prefitParamSlope_[iMomBin]->setConstant(slope_isConstant);
      prefitParamOffset_[iMomBin]->setConstant(offset_isConstant);
      prefitParamC_[iMomBin]->setConstant(C_isConstant);
      prefitParamMP2_[iMomBin]->setConstant(mp2_isConstant);
      prefitParamWidth2_[iMomBin]->setConstant(width2_isConstant);
      prefitParamMP3_[iMomBin]->setConstant(mp3_isConstant);
      prefitParamWidth3_[iMomBin]->setConstant(width3_isConstant);
      prefitParamX0_[iMomBin]->setConstant(x0_isConstant);
      prefitParamDeltaX1_[iMomBin]->setConstant(dx1_isConstant);
      prefitParamDeltaX2_[iMomBin]->setConstant(dx2_isConstant);

      std::cout << "fitParameter initial values:" << std::endl;
      std::cout << " prefitParamMP1      = " << prefitParamMP1_[iMomBin]->getVal()     << "," 
		<< " isConstant = " << prefitParamMP1_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamWidth1   = " << prefitParamWidth1_[iMomBin]->getVal()  << ","
		<< " isConstant = " << prefitParamWidth1_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamGMean    = " << prefitParamGMean_[iMomBin]->getVal()   << ","
		<< " isConstant = " << prefitParamGMean_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamGSigma   = " << prefitParamGSigma_[iMomBin]->getVal()  << ","
		<< " isConstant = " << prefitParamGSigma_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamAlpha    = " << prefitParamAlpha_[iMomBin]->getVal()   << ","
		<< " isConstant = " << prefitParamAlpha_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamSlope    = " << prefitParamSlope_[iMomBin]->getVal()   << ","
		<< " isConstant = " << prefitParamSlope_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamOffset   = " << prefitParamOffset_[iMomBin]->getVal()  << ","
		<< " isConstant = " << prefitParamOffset_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamC        = " << prefitParamC_[iMomBin]->getVal()       << ","
		<< " isConstant = " << prefitParamC_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamMP2      = " << prefitParamMP2_[iMomBin]->getVal()     << ","
		<< " isConstant = " << prefitParamMP2_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamWidth2   = " << prefitParamWidth2_[iMomBin]->getVal()  << ","
		<< " isConstant = " << prefitParamWidth2_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamMP3      = " << prefitParamMP3_[iMomBin]->getVal()     << ","
		<< " isConstant = " << prefitParamMP3_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamWidth3   = " << prefitParamWidth3_[iMomBin]->getVal()  << ","
		<< " isConstant = " << prefitParamWidth3_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamX0       = " << prefitParamX0_[iMomBin]->getVal()      << ","
		<< " isConstant = " << prefitParamX0_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamDeltaX1  = " << prefitParamDeltaX1_[iMomBin]->getVal() << ","
		<< " isConstant = " << prefitParamDeltaX1_[iMomBin]->isConstant() << std::endl;
      std::cout << " prefitParamDeltaX2  = " << prefitParamDeltaX2_[iMomBin]->getVal() << ","
		<< " isConstant = " << prefitParamDeltaX2_[iMomBin]->isConstant() << std::endl;

      TString modelName = Form("pdf_%s_%s_%s_%s", decayMode_string.Data(), momBinName.Data(), sep_->GetName(), label_.Data());
      prefitModel_[iMomBin] = 
	new TauDecayKinePdf2(modelName.Data(), modelName.Data(), 
			     *sepTimesMom_prefit, 
			     *prefitParamMP1_[iMomBin], *prefitParamWidth1_[iMomBin], 
			     *prefitParamGMean_[iMomBin], *prefitParamGSigma_[iMomBin], *prefitParamAlpha_[iMomBin], 
			     *prefitParamSlope_[iMomBin], *prefitParamOffset_[iMomBin], *prefitParamC_[iMomBin],
			     *prefitParamMP2_[iMomBin], *prefitParamWidth2_[iMomBin], 
			     *prefitParamMP3_[iMomBin], *prefitParamWidth3_[iMomBin], 
			     *prefitParamX0_[iMomBin], *prefitParamDeltaX1_[iMomBin], *prefitParamDeltaX2_[iMomBin]);

      prefitModel_[iMomBin]->disableAnalyticIntegration();

      TObjArray prefitParamConstraints;

      if ( !gmean_isConstant ) {
	RooConstVar* gmeanConstraint_value = 
	  new RooConstVar("gmeanConstraint_value", "gmeanConstraint_value", gmeanInitial);
	RooConstVar* gmeanConstraint_sigma =
	  new RooConstVar("gmeanConstraint_sigma", "gmeanConstraint_sigma", 0.25*gmeanInitial);
	RooGaussian* gmeanConstraint_pdf =
	  new RooGaussian("gmeanConstraint_pdf", "gmeanConstraint_pdf",
			  *prefitParamGMean_[iMomBin], *gmeanConstraint_value, *gmeanConstraint_sigma);
	prefitParamConstraints.Add(gmeanConstraint_pdf);
      }

      if ( !slope_isConstant ) {
	RooConstVar* slopeConstraint_value = 
	  new RooConstVar("slopeConstraint_value", "slopeConstraint_value", slopeInitial);
	RooConstVar* slopeConstraint_sigma =
	  new RooConstVar("slopeConstraint_sigma", "slopeConstraint_sigma", 10.*slopeInitial);
	RooGaussian* slopeConstraint_pdf =
	  new RooGaussian("slopeConstraint_pdf", "slopeConstraint_pdf",
			  *prefitParamSlope_[iMomBin], *slopeConstraint_value, *slopeConstraint_sigma);
	prefitParamConstraints.Add(slopeConstraint_pdf);
      }

      if ( !dx1_isConstant ) {
	RooConstVar* dx1Constraint_value =
	  new RooConstVar("dx1Constraint_value", "dx1Constraint_value", dx1Initial);
	RooConstVar* dx1Constraint_sigma =
	  new RooConstVar("dx1Constraint_sigma", "dx1Constraint_sigma", 0.5*dx1Initial);
	RooGaussian* dx1Constraint_pdf =
	  new RooGaussian("dx1Constraint_pdf", "dx1Constraint_pdf",
			  *prefitParamDeltaX1_[iMomBin], *dx1Constraint_value, *dx1Constraint_sigma);
	prefitParamConstraints.Add(dx1Constraint_pdf);
      }

      if ( decayMode_ == kElectron_Muon  ||
	   decayMode_ == kOneProngGt0Pi0 ||
	   decayMode_ == kThreeProng0Pi0 ) {
	RooConstVar* dx2Constraint_value =
	  new RooConstVar("dx2Constraint_value", "dx2Constraint_value", 25.);
	RooConstVar* dx2Constraint_sigma =
	  new RooConstVar("dx2Constraint_sigma", "dx2Constraint_sigma", 25.);
	RooGaussian* dx2Constraint_pdf =
	  new RooGaussian("dx2Constraint_pdf", "dx2Constraint_pdf",
			  *prefitParamDeltaX1_[iMomBin], *dx2Constraint_value, *dx2Constraint_sigma);
	prefitParamConstraints.Add(dx2Constraint_pdf);
      }

      RooLinkedList options;
      options.Add(new RooCmdArg(RooFit::Save(true)));
      //options.Add(new RooCmdArg(RooFit::SumW2Error(true)));
      options.Add(new RooCmdArg(RooFit::PrintLevel(-1)));
      options.Add(new RooCmdArg(RooFit::PrintEvalErrors(-1)));
      options.Add(new RooCmdArg(RooFit::Warnings(-1)));
      options.Add(new RooCmdArg(RooFit::ExternalConstraints(RooArgSet(prefitParamConstraints))));

//--- perform stand-alone fit of first Landau distribution
      if ( !(mp1_isConstant && width1_isConstant) ) {
	std::cout << "--> fitting Landau1 distribution..." << std::endl;

//--- set bin-contents of empty bins to one,
//    and include them in the fit
	int prePeak_bin_left = histogram->FindBin(prePeak_x_left);
	for ( int iBin = (prePeak_bin_left - 3); iBin <= prePeak_bin_left; ++iBin ) {
	  if ( iBin >= 1 && histogram->GetBinContent(iBin) < 1. ) {
	    histogram->SetBinContent(iBin, 1.);
	    histogram->SetBinError(iBin, 1.);
	  }
	}

	double sepTimesMomMin_landau1 = prePeak_x_left - 2.*TMath::Abs(prePeak_x_right - prePeak_x_left);
	double sepTimesMomMax_landau1 = prePeak_x_right;
	sepTimesMom_prefit->setRange("landau1", sepTimesMomMin_landau1, sepTimesMomMax_landau1);
	std::cout << " range = " << sepTimesMomMin_landau1 << ".." << sepTimesMomMax_landau1 << std::endl;

	RooLinkedList options_landau1(options);
	options_landau1.Add(new RooCmdArg(RooFit::Range("landau1")));

	if ( !mp1_isConstant    ) prefitParamMP1_[iMomBin]->setConstant(false);
	if ( !width1_isConstant ) prefitParamWidth1_[iMomBin]->setConstant(false);
	prefitParamGMean_[iMomBin]->setConstant(true);
	prefitParamGSigma_[iMomBin]->setConstant(true);
	prefitParamAlpha_[iMomBin]->setConstant(true);
	prefitParamSlope_[iMomBin]->setConstant(true);
        prefitParamOffset_[iMomBin]->setConstant(true);
	prefitParamC_[iMomBin]->setConstant(true);
	prefitParamMP2_[iMomBin]->setConstant(true);
	prefitParamWidth2_[iMomBin]->setConstant(true);
	prefitParamMP3_[iMomBin]->setConstant(true);
	prefitParamWidth3_[iMomBin]->setConstant(true);
	prefitParamX0_[iMomBin]->setConstant(true);
	prefitParamDeltaX1_[iMomBin]->setConstant(true);
	prefitParamDeltaX2_[iMomBin]->setConstant(true);

	RooFitResult* prefitResult_landau1 = prefitModel_[iMomBin]->fitTo(*datahist, options_landau1);
	delete prefitResult_landau1;

	printAllPrefitParameter(iMomBin);
      } else {
	std::cout << "--> skipping fit of Landau1 distribution..." << std::endl;
      }

//--- perform stand-alone fit of Gaussian distribution
      if ( !(gmean_isConstant && gsigma_isConstant && alpha_isConstant && slope_isConstant && offset_isConstant && C_isConstant) ) {
	std::cout << "--> fitting Gaussian distribution..." << std::endl;

	double sepTimesMomMin_gaussian = x0Initial;
	double sepTimesMomMax_gaussian = x0Initial + dx1Initial;
	sepTimesMom_prefit->setRange("gaussian", sepTimesMomMin_gaussian, sepTimesMomMax_gaussian);
	std::cout << " range = " << sepTimesMomMin_gaussian << ".." << sepTimesMomMax_gaussian << std::endl;

	RooLinkedList options_gaussian(options);
	options_gaussian.Add(new RooCmdArg(RooFit::Range("gaussian")));

	prefitParamMP1_[iMomBin]->setConstant(true);
	prefitParamWidth1_[iMomBin]->setConstant(true);
	if ( !gmean_isConstant  ) prefitParamGMean_[iMomBin]->setConstant(false);
	if ( !gsigma_isConstant ) prefitParamGSigma_[iMomBin]->setConstant(false);
	if ( !alpha_isConstant  ) prefitParamAlpha_[iMomBin]->setConstant(false);
	if ( !slope_isConstant  ) prefitParamSlope_[iMomBin]->setConstant(false);
	if ( !offset_isConstant ) prefitParamOffset_[iMomBin]->setConstant(false);
	if ( !C_isConstant      ) prefitParamC_[iMomBin]->setConstant(false);
	prefitParamMP2_[iMomBin]->setConstant(true);
	prefitParamWidth2_[iMomBin]->setConstant(true);
	prefitParamMP3_[iMomBin]->setConstant(true);
	prefitParamWidth3_[iMomBin]->setConstant(true);
	prefitParamX0_[iMomBin]->setConstant(true);
	prefitParamDeltaX1_[iMomBin]->setConstant(true);
	prefitParamDeltaX2_[iMomBin]->setConstant(true);

	RooFitResult* prefitResult_gaussian = prefitModel_[iMomBin]->fitTo(*datahist, options_gaussian);
	delete prefitResult_gaussian;

	printAllPrefitParameter(iMomBin);

	// CV: adjust x1 position to match Gaussian mean,
	//     relax condition on MP2 parameter
	if ( decayMode_ == kOneProng0Pi0 && isSelected_ ) {
	  //gmean_isConstant = true;
	  //slope_isConstant = true;
	  //offset_isConstant = true;
	  if ( prefitParamC_[iMomBin]->getVal() > epsilon ) {
            dx1Initial = TMath::Max(2.5*histogramMax_binWidth, prefitParamGMean_[iMomBin]->getVal() - prefitParamX0_[iMomBin]->getVal());
	    mp2Initial = prefitParamGMean_[iMomBin]->getVal();	    	    
	    mp2Min     = mp2Initial - 1.;
	    mp2Max     = mp2Initial + 0.25;
	    std::cout << "--> adjusting mp2Initial = " << mp2Initial << "," 
		      << " range = " << mp2Min << ".." << mp2Max << std::endl;
	    setRealVar_Value_Range(prefitParamMP2_[iMomBin], mp2Initial, mp2Min, mp2Max);
	  } else {
	    dx1Initial = TMath::Max(2.*histogramMax_binWidth, histogramMax_x - prefitParamX0_[iMomBin]->getVal());
	    //if ( dx1Initial < 2.5*histogramMax_binWidth ) dx1_isConstant = true;
	  }
	  dx1Min = TMath::Max(1.5*histogramMax_binWidth, dx1Initial - 0.25);
	  dx1Max = dx1Initial + 1.;
	  std::cout << "--> adjusting dx1Initial = " << dx1Initial << "," 
		    << " range = " << dx1Min << ".." << dx1Max << std::endl;
	  setRealVar_Value_Range(prefitParamDeltaX1_[iMomBin], dx1Initial, dx1Min, dx1Max);
	  if ( dx1_isConstant ) prefitParamDeltaX1_[iMomBin]->setConstant(true);
	}
      } else {
	std::cout << "--> skipping fit of Gaussian distribution..." << std::endl;
      }

      double x2 = x0Initial + dx1Initial + dx2Initial;
      if ( x2 > 25. ) x2 = 25.;

//--- perform stand-alone fit of first Landau distribution
      if ( !(mp2_isConstant && width2_isConstant) ) {
	std::cout << "--> fitting Landau2 distribution..." << std::endl;
	
	int iBin0 = histogram->FindBin(x0Initial + dx1Initial);
	const int numBinsRef = 5;
	int iBinRef = iBin0 + numBinsRef;
	std::cout << "iBinRef = " << iBinRef << ": position = " << histogram->GetBinCenter(iBinRef) << std::endl;

	double sepTimesMomMin_landau2 = histogram->GetBinCenter(iBin0);
	double sepTimesMomMax_landau2 = histogram->GetBinCenter(iBinRef);
	sepTimesMom_prefit->setRange("landau2", sepTimesMomMin_landau2, sepTimesMomMax_landau2);
	std::cout << " range = " << sepTimesMomMin_landau2 << ".." << sepTimesMomMax_landau2 << std::endl;

	RooLinkedList options_landau2(options);
	options_landau2.Add(new RooCmdArg(RooFit::Range("landau2")));

	prefitParamMP1_[iMomBin]->setConstant(true);
	prefitParamWidth1_[iMomBin]->setConstant(true);
	prefitParamGMean_[iMomBin]->setConstant(true);
	prefitParamGSigma_[iMomBin]->setConstant(true);
	prefitParamAlpha_[iMomBin]->setConstant(true);
	prefitParamSlope_[iMomBin]->setConstant(true);
	prefitParamOffset_[iMomBin]->setConstant(true);
	prefitParamC_[iMomBin]->setConstant(true);
	if ( !mp2_isConstant    ) prefitParamMP2_[iMomBin]->setConstant(false);
	if ( !width2_isConstant ) prefitParamWidth2_[iMomBin]->setConstant(false);
	prefitParamMP3_[iMomBin]->setConstant(true);
	prefitParamWidth3_[iMomBin]->setConstant(true);
	prefitParamX0_[iMomBin]->setConstant(true);
	prefitParamDeltaX1_[iMomBin]->setConstant(true);
	prefitParamDeltaX2_[iMomBin]->setConstant(true);

	RooFitResult* prefitResult_landau2 = prefitModel_[iMomBin]->fitTo(*datahist, options_landau2);
	delete prefitResult_landau2;
      
	if ( !dx2_isConstant ) {
	  // CV: return value of RooFitResult::minNll is unreliable,
	  //     because it may change by large amount in case Hesse matrix has negative diagonal elements
	  //    --> call dedicated function in order to compute FCN value a-la Minuit
	  double refFCN = compFCN(histogram, iBin0, iBinRef, 
				  prefitParamMP2_[iMomBin]->getVal(), prefitParamWidth2_[iMomBin]->getVal());
	  //std::cout << "refFCN = " << refFCN << std::endl;
	  double refFCNperDoF = refFCN/((iBinRef - iBin0) + 1 - 2);
	  //std::cout << "refFCNperDoF = " << refFCNperDoF << std::endl;
	  int iBin1 = histogram->GetNbinsX();
	  //std::cout << "iBin1 = " << iBin1 << ": position = " << histogram->GetBinCenter(iBin1) << std::endl;
	  int dBin = TMath::CeilNint(0.5*(histogram->GetNbinsX() - iBin0));
	  //std::cout << "dBin = " << dBin << std::endl;
	  double FCNperDoF = 0.;
	  do {	
	    sepTimesMom_prefit->setRange("landau2", histogram->GetBinCenter(iBin0), histogram->GetBinCenter(iBin1));
	    RooFitResult* prefitResult_landau2 = prefitModel_[iMomBin]->fitTo(*datahist, options_landau2);
	    delete prefitResult_landau2;	
	    double FCN = compFCN(histogram, iBin0, iBin1, 
				 prefitParamMP2_[iMomBin]->getVal(), prefitParamWidth2_[iMomBin]->getVal());
	    //std::cout << "FCN = " << FCN << std::endl;
	    FCNperDoF = FCN/((iBin1 - iBin0) + 1 - 2);
	    //std::cout << "FCNperDoF = " << FCNperDoF << std::endl;
	    if ( FCNperDoF < refFCNperDoF ) refFCNperDoF = FCNperDoF;
	    //if ( FCNperDoF > (1.25*refFCNperDoF) ) iBin1 = TMath::Max(iBinRef, iBin1 - dBin);
	    if ( FCNperDoF > (2.*refFCNperDoF) ) iBin1 = TMath::Max(iBinRef, iBin1 - dBin);
	    else iBin1 = TMath::Min(histogram->GetNbinsX(), iBin1 + dBin);
	    //std::cout << "iBin1 = " << iBin1 << ": position = " << histogram->GetBinCenter(iBin1) << std::endl;
	    dBin = TMath::CeilNint(0.5*dBin);
	    //std::cout << "dBin = " << dBin << std::endl;
	  } while ( dBin > 1 );
	  
	  //if ( FCNperDoF > (1.25*refFCNperDoF) ) iBin1 = TMath::Max(iBinRef, iBin1 - dBin);
	  if ( FCNperDoF > (2.*refFCNperDoF) ) iBin1 = TMath::Max(iBinRef, iBin1 - dBin);
	  iBin1 -= 3; // CV: "phenomenological" correction...
	  
	  double x2 = histogram->GetBinCenter(iBin1);
	  if ( iBin1 > (histogram->GetNbinsX() - 5) ) x2 = 25.;
	  std::cout << "x2 = " << x2 << std::endl;
	  
	  double dx2Initial = x2 - (x0Initial + dx1Initial);
	  double dx2Min     = dx2Initial - 2.*histogram->GetBinWidth(iBin1);
	  double dx2Max     = dx2Initial + 2.*histogram->GetBinWidth(iBin1);
	  if ( dx2Max > 25. ) dx2Max = 25.;
	  setRealVar_Value_Range(prefitParamDeltaX2_[iMomBin], dx2Initial, dx2Min, dx2Max);
	  setRealVar_Value_Range(prefitParamMP3_[iMomBin], 0.8*x2, 0.6*x2, x2);
	}

	printAllPrefitParameter(iMomBin);
      } else {
	std::cout << "--> skipping fit of Landau2 distribution..." << std::endl;
      }

//--- perform stand-alone fit of second Landau3 distribution
      if ( x2 < 25. && !(mp3_isConstant && width3_isConstant)) {
	std::cout << "--> fitting Landau3 distribution..." << std::endl;
	
	double sepTimesMomMin_landau3 = x2;
	double sepTimesMomMax_landau3 = 25.;
	sepTimesMom_prefit->setRange("landau3", sepTimesMomMin_landau3, sepTimesMomMax_landau3);
	std::cout << " range = " << sepTimesMomMin_landau3 << ".." << sepTimesMomMax_landau3 << std::endl;

	RooLinkedList options_landau3(options);
	options_landau3.Add(new RooCmdArg(RooFit::Range("landau3")));

	prefitParamMP1_[iMomBin]->setConstant(true);
	prefitParamWidth1_[iMomBin]->setConstant(true);
	prefitParamGMean_[iMomBin]->setConstant(true);
	prefitParamGSigma_[iMomBin]->setConstant(true);
	prefitParamAlpha_[iMomBin]->setConstant(true);	
	prefitParamSlope_[iMomBin]->setConstant(true);
	prefitParamOffset_[iMomBin]->setConstant(true);
	prefitParamC_[iMomBin]->setConstant(true);
	prefitParamMP2_[iMomBin]->setConstant(true);
	prefitParamWidth2_[iMomBin]->setConstant(true);
	if ( !mp3_isConstant    ) prefitParamMP3_[iMomBin]->setConstant(false);
	if ( !width3_isConstant ) prefitParamWidth3_[iMomBin]->setConstant(false);
	prefitParamX0_[iMomBin]->setConstant(true);
	prefitParamDeltaX1_[iMomBin]->setConstant(true);
	prefitParamDeltaX2_[iMomBin]->setConstant(true);
	
	RooFitResult* prefitResult_landau3 = prefitModel_[iMomBin]->fitTo(*datahist, options_landau3);
	delete prefitResult_landau3;
	
	printAllPrefitParameter(iMomBin);
      } else {
	std::cout << "--> skipping fit of Landau3 distribution..." << std::endl;
      }

//--- start combined fit of Gaussian + Landau + Landau2 model
      std::cout << "--> starting combined fit of Gaussian + Landau1 + Landau2 model..." << std::endl;
      
      double sepTimesMomMin_combined =  0.;
      double sepTimesMomMax_combined = 25.;
      //if ( decayMode_ == kOneProng0Pi0 ) {
      //  sepTimesMomMin_combined = prePeak_x_left - 2.*TMath::Abs(prePeak_x_right - prePeak_x_left);
      //}
      sepTimesMom_prefit->setRange("combined", sepTimesMomMin_combined, sepTimesMomMax_combined);

      RooLinkedList options_combined(options);
      options_combined.Add(new RooCmdArg(RooFit::Range("combined")));

      prefitParamMP1_[iMomBin]->setConstant(true);
      prefitParamWidth1_[iMomBin]->setConstant(true);
      if ( !gmean_isConstant  ) prefitParamGMean_[iMomBin]->setConstant(false);
      if ( !gsigma_isConstant ) prefitParamGSigma_[iMomBin]->setConstant(false);
      if ( !alpha_isConstant  ) prefitParamAlpha_[iMomBin]->setConstant(false);
      if ( !slope_isConstant  ) prefitParamSlope_[iMomBin]->setConstant(false);
      if ( !offset_isConstant ) prefitParamOffset_[iMomBin]->setConstant(false);
      if ( !C_isConstant      ) prefitParamC_[iMomBin]->setConstant(false);
      if ( !mp2_isConstant    ) prefitParamMP2_[iMomBin]->setConstant(false);
      if ( !width2_isConstant ) prefitParamWidth2_[iMomBin]->setConstant(false);
      if ( !mp3_isConstant    ) prefitParamMP3_[iMomBin]->setConstant(false);
      if ( !width3_isConstant ) prefitParamWidth3_[iMomBin]->setConstant(false);
      prefitParamX0_[iMomBin]->setConstant(true);
      if ( !dx1_isConstant    ) prefitParamDeltaX1_[iMomBin]->setConstant(false);
      if ( !dx2_isConstant    ) prefitParamDeltaX2_[iMomBin]->setConstant(false);

      RooFitResult* prefitResult = prefitModel_[iMomBin]->fitTo(*datahist, options_combined);
      std::cout << " prefit status = " << prefitResult->status() << " (converged = 0)" << std::endl;
      delete prefitResult;

//--- CV: map first Landau --> second Landau in case only one Landau distribution is needed to fit the TAUOLA prediction,
//        in order to make dx1 linearly increasing at low pt/energy
//       (TGraph fits will not converge otherwise)
      if ( (prefitParamX0_[iMomBin]->getVal() + prefitParamDeltaX1_[iMomBin]->getVal() + prefitParamDeltaX2_[iMomBin]->getVal()) > 24.5 &&
	   ((0.5*(momMin + momMax) < 75. && decayMode_string == "Electron_Muon")  || 
	    (0.5*(momMin + momMax) < 75. && decayMode_string == "ThreeProng0Pi0")) ) {
	std::cout << "mapping Landau1 --> Landau2." << std::endl;
	setRealVar_Value_Range(prefitParamMP3_[iMomBin], prefitParamMP2_[iMomBin]->getVal(), 
			       prefitParamMP2_[iMomBin]->getVal() - epsilon, prefitParamMP2_[iMomBin]->getVal() + epsilon);
	setRealVar_Value_Range(prefitParamWidth3_[iMomBin], prefitParamWidth2_[iMomBin]->getVal(), 
			       prefitParamWidth2_[iMomBin]->getVal() - epsilon, prefitParamWidth2_[iMomBin]->getVal() + epsilon);
	setRealVar_Value_Range(prefitParamDeltaX2_[iMomBin], 0., -epsilon, +epsilon);
      }

      printAllPrefitParameter(iMomBin);

      prefitModel_[iMomBin]->enableAnalyticIntegration();
      
      if ( prefitParamX0_[iMomBin]->getVal() > 0. ) {
	storePrefitResults(prefitParamMP1_[iMomBin], momMin, momMax, prefitResultMP1_);
	storePrefitResults(prefitParamWidth1_[iMomBin], momMin, momMax, prefitResultWidth1_);
      }
      storePrefitResults(prefitParamGMean_[iMomBin], momMin, momMax, prefitResultGMean_);
      storePrefitResults(prefitParamGSigma_[iMomBin], momMin, momMax, prefitResultGSigma_);
      storePrefitResults(prefitParamAlpha_[iMomBin], momMin, momMax, prefitResultAlpha_);
      storePrefitResults(prefitParamSlope_[iMomBin], momMin, momMax, prefitResultSlope_);
      storePrefitResults(prefitParamOffset_[iMomBin], momMin, momMax, prefitResultOffset_);
      storePrefitResults(prefitParamC_[iMomBin], momMin, momMax, prefitResultC_);
      storePrefitResults(prefitParamMP2_[iMomBin], momMin, momMax, prefitResultMP2_);
      storePrefitResults(prefitParamWidth2_[iMomBin], momMin, momMax, prefitResultWidth2_);
      if ( (prefitParamX0_[iMomBin]->getVal() 
           + prefitParamDeltaX1_[iMomBin]->getVal() + prefitParamDeltaX2_[iMomBin]->getVal()) < 24.5 ) {
	storePrefitResults(prefitParamMP3_[iMomBin], momMin, momMax, prefitResultMP3_);
	storePrefitResults(prefitParamWidth3_[iMomBin], momMin, momMax, prefitResultWidth3_);
      }
      storePrefitResults(prefitParamX0_[iMomBin], momMin, momMax, prefitResultX0_);
      storePrefitResults(prefitParamDeltaX1_[iMomBin], momMin, momMax, prefitResultDeltaX1_);
      storePrefitResults(prefitParamDeltaX2_[iMomBin], momMin, momMax, prefitResultDeltaX2_);

      for ( unsigned iScale = 0; iScale < 2; ++iScale ) {
	canvas->Clear();

	if ( iScale == 0 ) canvas->SetLogy(true);
	else               canvas->SetLogy(false);
    
	TString frameTitle = Form("%s %s: P_{T} = %2.0f..%2.0f GeV", sep_->GetName(), decayMode_string.Data(), momMin, momMax);
	RooPlot* frame = sepTimesMom_prefit->frame(RooFit::Title(frameTitle.Data()), RooFit::Bins(250));

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
	if ( iScale == 0 ) outputFileName_i.ReplaceAll(".", "_log.");
	else               outputFileName_i.ReplaceAll(".", "_linear.");
	
	canvas->Update();
	
	canvas->SaveAs(outputFileName_i.Data());
	canvas->SaveAs(outputFileName_i.ReplaceAll(".eps", ".root").Data());
      }

      delete datahist;
    }

    delete canvas;

    delete inputFile;
   
    delete sepTimesMom_prefit;

    fitTF1(momParamMP1_,     setMinMP1_,     setMaxMP1_,     prefitResultMP1_,     prefitCoeffValMP1_,     prefitCoeffErrMP1_);
    fitTF1(momParamWidth1_,  setMinWidth1_,  setMaxWidth1_,  prefitResultWidth1_,  prefitCoeffValWidth1_,  prefitCoeffErrWidth1_);
    fitTF1(momParamGMean_,   setMinGMean_,   setMaxGMean_,   prefitResultGMean_,   prefitCoeffValGMean_,   prefitCoeffErrGMean_);
    fitTF1(momParamGSigma_,  setMinGSigma_,  setMaxGSigma_,  prefitResultGSigma_,  prefitCoeffValGSigma_,  prefitCoeffErrGSigma_);
    fitTF1(momParamAlpha_,   setMinAlpha_,   setMaxAlpha_,   prefitResultAlpha_,   prefitCoeffValAlpha_,   prefitCoeffErrAlpha_);
    fitTF1(momParamSlope_,   setMinSlope_,   setMaxSlope_,   prefitResultSlope_,   prefitCoeffValSlope_,   prefitCoeffErrSlope_);
    fitTF1(momParamOffset_,  setMinOffset_,  setMaxOffset_,  prefitResultOffset_,  prefitCoeffValOffset_,  prefitCoeffErrOffset_);
    fitTF1(momParamC_,       setMinC_,       setMaxC_,       prefitResultC_,       prefitCoeffValC_,       prefitCoeffErrC_);
    fitTF1(momParamMP2_,     setMinMP2_,     setMaxMP2_,     prefitResultMP2_,     prefitCoeffValMP2_,     prefitCoeffErrMP2_);
    fitTF1(momParamWidth2_,  setMinWidth2_,  setMaxWidth2_,  prefitResultWidth2_,  prefitCoeffValWidth2_,  prefitCoeffErrWidth2_);
    fitTF1(momParamMP3_,     setMinMP3_,     setMaxMP3_,     prefitResultMP3_,     prefitCoeffValMP3_,     prefitCoeffErrMP3_);
    fitTF1(momParamWidth3_,  setMinWidth3_,  setMaxWidth3_,  prefitResultWidth3_,  prefitCoeffValWidth3_,  prefitCoeffErrWidth3_);
    fitTF1(momParamX0_,      setMinX0_,      setMaxX0_,      prefitResultX0_,      prefitCoeffValX0_,      prefitCoeffErrX0_);
    fitTF1(momParamDeltaX1_, setMinDeltaX1_, setMaxDeltaX1_, prefitResultDeltaX1_, prefitCoeffValDeltaX1_, prefitCoeffErrDeltaX1_);
    fitTF1(momParamDeltaX2_, setMinDeltaX2_, setMaxDeltaX2_, prefitResultDeltaX2_, prefitCoeffValDeltaX2_, prefitCoeffErrDeltaX2_);
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

    fitParamMP1_     = buildMomDependentFitParameter("mp1", fitParamCoeffMP1_, 
						     prefitCoeffValMP1_, prefitCoeffErrMP1_, momParamMP1_->GetTitle());
    fitParamWidth1_  = buildMomDependentFitParameter("width1", fitParamCoeffWidth1_, 
						     prefitCoeffValWidth1_, prefitCoeffErrWidth1_, momParamWidth1_->GetTitle());
    fitParamGMean_   = buildMomDependentFitParameter("gmean", fitParamCoeffGMean_, 
						     prefitCoeffValGMean_, prefitCoeffErrGMean_, momParamGMean_->GetTitle());
    fitParamGSigma_  = buildMomDependentFitParameter("gsigma", fitParamCoeffGSigma_, 
						     prefitCoeffValGSigma_, prefitCoeffErrGSigma_, momParamGSigma_->GetTitle());
    fitParamAlpha_   = buildMomDependentFitParameter("alpha", fitParamCoeffAlpha_, 
						     prefitCoeffValAlpha_, prefitCoeffErrAlpha_, momParamAlpha_->GetTitle());
    fitParamSlope_   = buildMomDependentFitParameter("slope", fitParamCoeffSlope_, 
						     prefitCoeffValSlope_, prefitCoeffErrSlope_, momParamSlope_->GetTitle());
    fitParamOffset_  = buildMomDependentFitParameter("offset", fitParamCoeffOffset_, 
						     prefitCoeffValOffset_, prefitCoeffErrOffset_, momParamOffset_->GetTitle());
    fitParamC_       = buildMomDependentFitParameter("C", fitParamCoeffC_, 
						     prefitCoeffValC_, prefitCoeffErrC_, momParamC_->GetTitle());
    fitParamMP2_     = buildMomDependentFitParameter("mp2", fitParamCoeffMP2_, 
						     prefitCoeffValMP2_, prefitCoeffErrMP2_, momParamMP2_->GetTitle());
    fitParamWidth2_  = buildMomDependentFitParameter("width2", fitParamCoeffWidth2_, 
						     prefitCoeffValWidth2_, prefitCoeffErrWidth2_, momParamWidth2_->GetTitle());
    fitParamMP3_     = buildMomDependentFitParameter("mp3", fitParamCoeffMP3_, 
						     prefitCoeffValMP3_, prefitCoeffErrMP3_, momParamMP3_->GetTitle());
    fitParamWidth3_  = buildMomDependentFitParameter("width3", fitParamCoeffWidth3_,
						     prefitCoeffValWidth3_, prefitCoeffErrWidth3_, momParamWidth3_->GetTitle());
    fitParamX0_      = buildMomDependentFitParameter("x0", fitParamCoeffX0_, 
						     prefitCoeffValX0_, prefitCoeffErrX0_, momParamX0_->GetTitle());
    fitParamDeltaX1_ = buildMomDependentFitParameter("dx1", fitParamCoeffDeltaX1_, 
						     prefitCoeffValDeltaX1_, prefitCoeffErrDeltaX1_, momParamDeltaX1_->GetTitle());
    fitParamDeltaX2_ = buildMomDependentFitParameter("dx2", fitParamCoeffDeltaX2_, 
						     prefitCoeffValDeltaX2_, prefitCoeffErrDeltaX2_, momParamDeltaX2_->GetTitle());
        
    TString modelName = Form("pdf_%s_%s_%s_%s", decayMode_string.Data(), "AllMom", sep_->GetName(), label_.Data());
    fitModel_ = 
      new TauDecayKinePdf2(modelName.Data(), modelName.Data(), 
			   *sepTimesMom_, 
			   *fitParamMP1_, *fitParamWidth1_, 
			   *fitParamGMean_, *fitParamGSigma_, *fitParamAlpha_, *fitParamSlope_, *fitParamOffset_, *fitParamC_,
			   *fitParamMP2_, *fitParamWidth2_, 
			   *fitParamMP3_, *fitParamWidth3_, 
			   *fitParamX0_, *fitParamDeltaX1_, *fitParamDeltaX2_);
  
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

    storeFitResults(fitParamCoeffMP1_, fitCoeffValMP1_, fitCoeffErrMP1_);
    storeFitResults(fitParamCoeffWidth1_, fitCoeffValWidth1_, fitCoeffErrWidth1_);
    storeFitResults(fitParamCoeffGMean_, fitCoeffValGMean_, fitCoeffErrGMean_);
    storeFitResults(fitParamCoeffGSigma_, fitCoeffValGSigma_, fitCoeffErrGSigma_);
    storeFitResults(fitParamCoeffAlpha_, fitCoeffValAlpha_, fitCoeffErrAlpha_);
    storeFitResults(fitParamCoeffSlope_, fitCoeffValSlope_, fitCoeffErrSlope_);
    storeFitResults(fitParamCoeffOffset_, fitCoeffValOffset_, fitCoeffErrOffset_);
    storeFitResults(fitParamCoeffC_, fitCoeffValC_, fitCoeffErrC_);
    storeFitResults(fitParamCoeffMP2_, fitCoeffValMP2_, fitCoeffErrMP2_);
    storeFitResults(fitParamCoeffWidth2_, fitCoeffValWidth2_, fitCoeffErrWidth2_);
    storeFitResults(fitParamCoeffMP3_, fitCoeffValMP3_, fitCoeffErrMP3_);
    storeFitResults(fitParamCoeffWidth3_, fitCoeffValWidth3_, fitCoeffErrWidth3_);
    storeFitResults(fitParamCoeffX0_, fitCoeffValX0_, fitCoeffErrX0_);
    storeFitResults(fitParamCoeffDeltaX1_, fitCoeffValDeltaX1_, fitCoeffErrDeltaX1_);
    storeFitResults(fitParamCoeffDeltaX2_, fitCoeffValDeltaX2_, fitCoeffErrDeltaX2_);

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
      RooPlot* frame = sepTimesMom_->frame(RooFit::Title(frameTitle.Data()), RooFit::Bins(250));
      
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

  bool isAngle_;
  bool isDeltaR_;

  bool isAll_;
  bool isSelected_;

  std::vector<RooRealVar*> prefitParamMP1_;
  std::vector<RooRealVar*> prefitParamWidth1_;
  std::vector<RooRealVar*> prefitParamGMean_;
  std::vector<RooRealVar*> prefitParamGSigma_;
  std::vector<RooRealVar*> prefitParamAlpha_;
  std::vector<RooRealVar*> prefitParamSlope_;
  std::vector<RooRealVar*> prefitParamOffset_;
  std::vector<RooRealVar*> prefitParamC_;
  std::vector<RooRealVar*> prefitParamMP2_;
  std::vector<RooRealVar*> prefitParamWidth2_;
  std::vector<RooRealVar*> prefitParamMP3_;
  std::vector<RooRealVar*> prefitParamWidth3_;
  std::vector<RooRealVar*> prefitParamX0_;
  std::vector<RooRealVar*> prefitParamDeltaX1_;
  std::vector<RooRealVar*> prefitParamDeltaX2_;
  std::vector<TauDecayKinePdf2*> prefitModel_;

  TGraphErrors* prefitResultMP1_;
  TGraphErrors* prefitResultWidth1_;
  TGraphErrors* prefitResultGMean_;
  TGraphErrors* prefitResultGSigma_;
  TGraphErrors* prefitResultAlpha_;
  TGraphErrors* prefitResultSlope_;
  TGraphErrors* prefitResultOffset_;
  TGraphErrors* prefitResultC_;
  TGraphErrors* prefitResultMP2_;
  TGraphErrors* prefitResultWidth2_;
  TGraphErrors* prefitResultMP3_;
  TGraphErrors* prefitResultWidth3_;
  TGraphErrors* prefitResultX0_;
  TGraphErrors* prefitResultDeltaX1_;
  TGraphErrors* prefitResultDeltaX2_;

  TString momParametrizationFormulaMP1_;
  bool setMinMP1_;
  bool setMaxMP1_;
  TString momParametrizationFormulaWidth1_;
  bool setMinWidth1_;
  bool setMaxWidth1_;
  TString momParametrizationFormulaGMean_;
  bool setMinGMean_;
  bool setMaxGMean_;
  TString momParametrizationFormulaGSigma_;
  bool setMinGSigma_;
  bool setMaxGSigma_;
  TString momParametrizationFormulaAlpha_;
  bool setMinAlpha_;
  bool setMaxAlpha_;
  TString momParametrizationFormulaSlope_;
  bool setMinSlope_;
  bool setMaxSlope_;
  TString momParametrizationFormulaOffset_;
  bool setMinOffset_;
  bool setMaxOffset_;
  TString momParametrizationFormulaC_;
  bool setMinC_;
  bool setMaxC_;
  TString momParametrizationFormulaMP2_;
  bool setMinMP2_;
  bool setMaxMP2_;  
  TString momParametrizationFormulaWidth2_;
  bool setMinWidth2_;
  bool setMaxWidth2_;
  TString momParametrizationFormulaMP3_;
  bool setMinMP3_;
  bool setMaxMP3_;
  TString momParametrizationFormulaWidth3_;
  bool setMinWidth3_;
  bool setMaxWidth3_;
  TString momParametrizationFormulaX0_;
  bool setMinX0_;
  bool setMaxX0_;
  TString momParametrizationFormulaDeltaX1_;
  bool setMinDeltaX1_;
  bool setMaxDeltaX1_;
  TString momParametrizationFormulaDeltaX2_;
  bool setMinDeltaX2_;
  bool setMaxDeltaX2_;  

  TF1* momParamMP1_;
  TF1* momParamWidth1_;
  TF1* momParamGMean_;
  TF1* momParamGSigma_;
  TF1* momParamAlpha_;
  TF1* momParamSlope_;
  TF1* momParamOffset_;
  TF1* momParamC_;
  TF1* momParamMP2_;
  TF1* momParamWidth2_;
  TF1* momParamMP3_;
  TF1* momParamWidth3_;
  TF1* momParamX0_;
  TF1* momParamDeltaX1_;
  TF1* momParamDeltaX2_;

  std::vector<double> prefitCoeffValMP1_;
  std::vector<double> prefitCoeffErrMP1_;
  std::vector<double> prefitCoeffValWidth1_;
  std::vector<double> prefitCoeffErrWidth1_;
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
  std::vector<double> prefitCoeffValMP2_;
  std::vector<double> prefitCoeffErrMP2_;
  std::vector<double> prefitCoeffValWidth2_;
  std::vector<double> prefitCoeffErrWidth2_;
  std::vector<double> prefitCoeffValMP3_;
  std::vector<double> prefitCoeffErrMP3_;
  std::vector<double> prefitCoeffValWidth3_;
  std::vector<double> prefitCoeffErrWidth3_;
  std::vector<double> prefitCoeffValX0_;
  std::vector<double> prefitCoeffErrX0_;
  std::vector<double> prefitCoeffValDeltaX1_;
  std::vector<double> prefitCoeffErrDeltaX1_;
  std::vector<double> prefitCoeffValDeltaX2_;
  std::vector<double> prefitCoeffErrDeltaX2_;
  
  std::vector<RooRealVar*> fitParamCoeffMP1_;
  RooAbsReal* fitParamMP1_;
  std::vector<RooRealVar*> fitParamCoeffWidth1_;
  RooAbsReal* fitParamWidth1_;
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
  std::vector<RooRealVar*> fitParamCoeffMP2_;
  RooAbsReal* fitParamMP2_;
  std::vector<RooRealVar*> fitParamCoeffWidth2_;
  RooAbsReal* fitParamWidth2_;
  std::vector<RooRealVar*> fitParamCoeffMP3_;
  RooAbsReal* fitParamMP3_;
  std::vector<RooRealVar*> fitParamCoeffWidth3_;
  RooAbsReal* fitParamWidth3_;
  std::vector<RooRealVar*> fitParamCoeffX0_;
  RooAbsReal* fitParamX0_;
  std::vector<RooRealVar*> fitParamCoeffDeltaX1_;
  RooAbsReal* fitParamDeltaX1_;
  std::vector<RooRealVar*> fitParamCoeffDeltaX2_;
  RooAbsReal* fitParamDeltaX2_;
  TauDecayKinePdf2* fitModel_;

  std::vector<double> fitCoeffValMP1_;
  std::vector<double> fitCoeffErrMP1_;
  std::vector<double> fitCoeffValWidth1_;
  std::vector<double> fitCoeffErrWidth1_;
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
  std::vector<double> fitCoeffValMP2_;
  std::vector<double> fitCoeffErrMP2_;
  std::vector<double> fitCoeffValWidth2_;
  std::vector<double> fitCoeffErrWidth2_;
  std::vector<double> fitCoeffValMP3_;
  std::vector<double> fitCoeffErrMP3_;
  std::vector<double> fitCoeffValWidth3_;
  std::vector<double> fitCoeffErrWidth3_;
  std::vector<double> fitCoeffValX0_;
  std::vector<double> fitCoeffErrX0_;
  std::vector<double> fitCoeffValDeltaX1_;
  std::vector<double> fitCoeffErrDeltaX1_;
  std::vector<double> fitCoeffValDeltaX2_;
  std::vector<double> fitCoeffErrDeltaX2_;
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
  RooRealVar leg1VisInvisAngleLabTimesEnergy("leg1VisInvisAngleLabTimesEnergy", "leg1VisInvisAngleLabTimesEnergy", 0., 25.);
  RooRealVar leg1VisInvisDeltaRLabTimesPt("leg1VisInvisDeltaRLabTimesPt", "leg1VisInvisDeltaRLabTimesPt", 0., 25.);

  leg1Energy.setBins(350);
  leg1Pt.setBins(250);
  leg1VisInvisAngleLab.setBins(250);
  leg1VisInvisDeltaRLab.setBins(250);
  leg1VisInvisAngleLabTimesEnergy.setBins(250);
  leg1VisInvisDeltaRLabTimesPt.setBins(250);

  RooRealVar leg2Energy("leg2Energy", "leg2Energy", 0., 350.);
  RooRealVar leg2Pt("leg2Pt", "leg2Pt", 0., 250.);
  RooRealVar leg2VisPt("leg2VisPt", "leg2VisPt", 0., 200.);
  RooRealVar leg2VisEta("leg2VisEta", "leg2VisEta", -2.5, +2.5);
  RooRealVar leg2DecayMode("leg2DecayMode", "leg2DecayMode", -1.5, +11.5);
  RooRealVar leg2VisInvisAngleLab("leg2VisInvisAngleLab", "leg2VisInvisAngleLab", 0., 1.);
  RooRealVar leg2VisInvisDeltaRLab("leg2VisInvisDeltaRLab", "leg2VisInvisDeltaRLab", 0., 1.);
  RooRealVar leg2VisInvisAngleLabTimesEnergy("leg2VisInvisAngleLabTimesEnergy", "leg2VisInvisAngleLabTimesEnergy", 0., 25.);
  RooRealVar leg2VisInvisDeltaRLabTimesPt("leg2VisInvisDeltaRLabTimesPt", "leg2VisInvisDeltaRLabTimesPt", 0., 25.);  

  leg2Energy.setBins(350);
  leg2Pt.setBins(250);
  leg2VisInvisAngleLab.setBins(250);
  leg2VisInvisDeltaRLab.setBins(250);
  leg2VisInvisAngleLabTimesEnergy.setBins(250);
  leg2VisInvisDeltaRLabTimesPt.setBins(250);

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
      
    TArrayD momBinning_all = getBinningMom(decayMode_string, "all");
    TArrayD momBinning_selected = getBinningMom(decayMode_string, "selected");
    TArrayD sepBinning = getBinningSepTimesMom(decayMode_string, selection_string);
    
    /*****************************************************************************
     ******** Run fit for leg1 with no cuts on visible decay products applied ****
     *****************************************************************************/

    if ( runLeg1 && runAll && runDeltaR ) {
      TString fitManagerName_leg1_dR_all = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg1", "dR", "all");
      if ( fitResults.find(fitManagerName_leg1_dR_all.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg1_dR_all.Data()] = 
	  new fitManager(&leg1Pt, &leg1VisInvisDeltaRLab, &leg1VisInvisDeltaRLabTimesPt, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg1", "dR", "all"), momBinning_all);
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
			 Form("%s_%s_%s", "leg1", "angle", "all"), momBinning_all);
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
			 Form("%s_%s_%s", "leg1", "dR", "selected1"), momBinning_selected);
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
			 Form("%s_%s_%s", "leg1", "angle", "selected1"), momBinning_selected);
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
			 Form("%s_%s_%s", "leg2", "dR", "all"), momBinning_all);
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
			 Form("%s_%s_%s", "leg2", "angle", "all"), momBinning_all);
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
			 Form("%s_%s_%s", "leg2", "dR", "selected1"), momBinning_selected);
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
			 Form("%s_%s_%s", "leg2", "angle", "selected1"), momBinning_selected);
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
