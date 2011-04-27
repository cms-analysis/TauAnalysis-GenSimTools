
#include "TauAnalysis/FittingTools/interface/SmoothLandau_x_GaussPdf.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooLandau.h"
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

enum { kElectron_Muon, 
       kOneProng0Pi0, kOneProng1Pi0, kOneProng2Pi0, kOneProngGt0Pi0, 
       kThreeProng0Pi0, kThreeProng1Pi0 };

TString getDecayMode_string(int decayMode)
{
  TString decayMode_string = "undefined";
  if      ( decayMode == kElectron_Muon  ) decayMode_string = "Electron_Muon";
  else if ( decayMode == kOneProng0Pi0   ) decayMode_string = "OneProng0Pi0";
  else if ( decayMode == kOneProng1Pi0   ) decayMode_string = "OneProng1Pi0"; 
  else if ( decayMode == kOneProng2Pi0   ) decayMode_string = "OneProng2Pi0"; 
  else if ( decayMode == kOneProngGt0Pi0 ) decayMode_string = "OneProngGt0Pi0";
  else if ( decayMode == kThreeProng0Pi0 ) decayMode_string = "ThreeProng0Pi0";
  else if ( decayMode == kThreeProng1Pi0 ) decayMode_string = "ThreeProng1Pi0";
  return decayMode_string;
}

double encodeDecayMode(int decayMode)
{
  double decayMode_encoded = -1000;
  if      ( decayMode == kElectron_Muon  ) decayMode_encoded = -1.;
  else if ( decayMode == kOneProng0Pi0   ) decayMode_encoded =  0.;
  else if ( decayMode == kOneProng1Pi0   ) decayMode_encoded =  1.;
  else if ( decayMode == kOneProng2Pi0   ) decayMode_encoded =  2.;
  else if ( decayMode == kThreeProng0Pi0 ) decayMode_encoded = 10.;
  else if ( decayMode == kThreeProng1Pi0 ) decayMode_encoded = 11.;
  return decayMode_encoded;
}

void pdfTimingTest(RooAbsPdf* pdf, RooRealVar* mom, RooRealVar* sep)
{
  TRandom3 rnd;

  TBenchmark benchmark;
  benchmark.Start("pdfTimingTest");

  double pdfValueSum = 0.;

  const unsigned numCalls = 10000;
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

void storePrefitResults(RooRealVar* fitParameter, double ptMin, double ptMax, TGraphErrors* graph)
{
  std::cout << "<storePrefitResults>:" << std::endl;

  int iPoint = graph->GetN();

  double pt = 0.5*(ptMin + ptMax);

  graph->SetPoint(iPoint, pt, fitParameter->getVal());

  graph->SetPointError(iPoint, 0.5*TMath::Abs(ptMax - ptMin), fitParameter->getError());
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
  fitManager(RooRealVar* mom, RooRealVar* dR, int decayMode, const TString& label,unsigned numPtBins, const double* ptBinning)
    : mom_(mom),
      dR_(dR),
      decayMode_(decayMode),
      label_(label),
      ptBinning_(numPtBins + 1, ptBinning),
      prefitParamLandauMP_(numPtBins),
      prefitParamLandauWidth_(numPtBins),
      prefitParamLandauScaleFactor_(numPtBins),
      prefitLandauScale_(numPtBins),
      prefitLandau_(numPtBins), 
      prefitParamGaussianMean_(numPtBins),
      prefitParamGaussianSigma_(numPtBins),
      prefitGaussian_(numPtBins), 
      prefitParamMix_(numPtBins),
      prefitModel_(numPtBins),
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
    double ptMin = ptBinning_[0];
    double ptMax = ptBinning_[ptBinning_.GetSize() - 1];
    TF1* tf1 = new TF1(tf1Name.Data(), formula, ptMin, ptMax);
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

  RooRealVar* buildPrefitParameter(const TString& name, const TString& ptBinName, double startValue, double min, double max)
  {
    std::cout << "<fitManager::buildPrefitParameter>:" << std::endl;

    TString decayMode_string = getDecayMode_string(decayMode_);
    TString fitParameterName = 
      Form("%s_%s_%s_%s_%s", name.Data(), ptBinName.Data(), dR_->GetName(), decayMode_string.Data(), label_.Data());
    RooRealVar* fitParameter = new RooRealVar(fitParameterName.Data(), fitParameterName.Data(), startValue, min, max);
    return fitParameter;
  }

  RooConstVar* buildConstPrefitParameter(const TString& name, const TString& ptBinName, double value)
  {
    std::cout << "<buildConstPrefitParameter>:" << std::endl;

    TString decayMode_string = getDecayMode_string(decayMode_);
    TString fitParameterName = 
      Form("%s_%s_%s_%s_%s", name.Data(), ptBinName.Data(), dR_->GetName(), decayMode_string.Data(), label_.Data());
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
  
    for ( int iPtBin = 0; iPtBin < (ptBinning_.GetSize() - 1); ++iPtBin ) {
      double ptMin = ptBinning_[iPtBin];
      double ptMax = ptBinning_[iPtBin + 1];
    
      TString histogramName = Form("%s_%sPt%2.0fto%2.0f", decayMode_string.Data(), label.Data(), ptMin, ptMax);      
      TH1* histogram = (TH1*)inputFile->Get(TString(inputDirName_full.Data()).Append("/").Append(histogramName.Data()));
      std::cout << " histogramName = " << histogramName.Data() << ": histogram = " << histogram << std::endl;

      TString datahistName = TString(histogramName).Append("_datahist");
      RooDataHist* datahist = new RooDataHist(datahistName.Data(), datahistName.Data(), RooArgList(*dR_), histogram);
      std::cout << " datahist = " << datahist << std::endl;

      TString ptBinName = Form("Pt%2.0fto%2.0f", ptMin, ptMax);

      double histogramMean = histogram->GetMean();
      if ( histogramMean > 0.15 ) histogramMean = 0.15;
      std::cout << " histogramMean = " << histogramMean << std::endl;
      double histogramRMS = histogram->GetRMS();
      if ( histogramRMS  > 0.15 ) histogramRMS  = 0.15;
      std::cout << " histogramRMS = " << histogramRMS << std::endl;
      
      prefitParamLandauMP_[iPtBin] = buildPrefitParameter("LandauMP", ptBinName, 1.2*histogramMean, 0., 0.50);
      prefitParamLandauWidth_[iPtBin] = buildPrefitParameter("LandauWidth", ptBinName, histogramRMS, 1.e-3, 0.50);
      prefitParamLandauScaleFactor_[iPtBin] = buildPrefitParameter("LandauScaleFactor", ptBinName, 1., 0.5, 2.0);
      TString LandauScaleName =
	Form("LandauScale_%s_%s_%s_%s", decayMode_string.Data(), ptBinName.Data(), dR_->GetName(), label_.Data());
      prefitLandauScale_[iPtBin] = 
	new RooFormulaVar(LandauScaleName.Data(), LandauScaleName.Data(), "@0/@1", RooArgSet(*dR_, *prefitParamLandauScaleFactor_[iPtBin]));
      prefitParamLandauLocation = buildPrefitParameter("LandauScaleLocation", ptBinName, 1., 0.5, 2.0);
      TString LandauName =
	Form("Landau_%s_%s_%s_%s", decayMode_string.Data(), ptBinName.Data(), dR_->GetName(), label_.Data());
      prefitLandau_[iPtBin] = 
	//	new RooLandau(LandauName.Data(), LandauName.Data(), *prefitLandauScale_[iPtBin],
	//	              *prefitParamLandauMP_[iPtBin], *prefitParamLandauWidth_[iPtBin]);
      new RooGenericPdf((LandauName.Data(), LandauName.Data(),
			 "TMath::Gaus(@0, @1, @2)*(1.0 + TMath::Erf(@3*((@0 - @1)/@2)/TMath::Sqrt(2)))",
			 RooArgList(*dR_, *prefitParamLandauMP_[iPtBin], *prefitParamLandauWidth_[iPtBin], prefitParamLandauLocation)));

*prefitParamSkewedGaussianMean_[iPtBin], *prefitParamSkewedGaussianSigma_[iPtBin], 
				    *prefitParamSkewedGaussianAlpha_[iPtBin])); 

      prefitParamGaussianMean_[iPtBin] = buildPrefitParameter("GaussianMean", ptBinName, 0.8*histogramMean, 0., 0.25);
      prefitParamGaussianSigma_[iPtBin] = buildPrefitParameter("GaussianSigma", ptBinName, histogramRMS, 1.e-3, 0.25);
      TString GaussianName = 
	Form("Gaussian_%s_%s_%s_%s", decayMode_string.Data(), ptBinName.Data(), dR_->GetName(), label_.Data());
      prefitGaussian_[iPtBin] = 
	new RooGaussian(GaussianName.Data(), GaussianName.Data(), *dR_,
      			*prefitParamGaussianMean_[iPtBin], *prefitParamGaussianSigma_[iPtBin]);
      
      prefitParamMix_[iPtBin] = buildPrefitParameter("mix", ptBinName, 0.50, 0., 1.);

      TString modelName = Form("pdf_%s_%s_%s_%s", decayMode_string.Data(), ptBinName.Data(), dR_->GetName(), label_.Data());
      prefitModel_[iPtBin] = 
	new RooAddPdf(modelName.Data(), modelName.Data(), 
		      *prefitLandau_[iPtBin], *prefitGaussian_[iPtBin], *prefitParamMix_[iPtBin]);

      RooLinkedList options;
      options.Add(new RooCmdArg(RooFit::Save(true)));
      options.Add(new RooCmdArg(RooFit::PrintLevel(-1)));
      options.Add(new RooCmdArg(RooFit::PrintEvalErrors(-1)));
      options.Add(new RooCmdArg(RooFit::Warnings(-1)));

      RooFitResult* prefitResult = prefitModel_[iPtBin]->fitTo(*datahist, options);
      std::cout << " prefit status = " << prefitResult->status() << " (converged = 0)" << std::endl;
      delete prefitResult;

      printPrefitParameter("LandauMP", prefitParamLandauMP_[iPtBin]);
      printPrefitParameter("LandauWidth", prefitParamLandauWidth_[iPtBin]);
      printPrefitParameter("LandauScaleFactor", prefitParamLandauScaleFactor_[iPtBin]);
      printPrefitParameter("GaussianMean", prefitParamGaussianMean_[iPtBin]);
      printPrefitParameter("GaussianSigma", prefitParamGaussianSigma_[iPtBin]);
      printPrefitParameter("mix", prefitParamMix_[iPtBin]);     

      storePrefitResults(prefitParamLandauMP_[iPtBin], ptMin, ptMax, prefitResultLandauMP_);
      storePrefitResults(prefitParamLandauWidth_[iPtBin], ptMin, ptMax, prefitResultLandauWidth_);
      storePrefitResults(prefitParamLandauScaleFactor_[iPtBin], ptMin, ptMax, prefitResultLandauScaleFactor_);
      storePrefitResults(prefitParamGaussianMean_[iPtBin], ptMin, ptMax, prefitResultGaussianMean_);
      storePrefitResults(prefitParamGaussianSigma_[iPtBin], ptMin, ptMax, prefitResultGaussianSigma_);
      storePrefitResults(prefitParamMix_[iPtBin], ptMin, ptMax, prefitResultMix_);
         
      canvas->Clear();
      canvas->SetLogy();
    
      TString frameTitle = Form("%s %s: P_{T} = %2.0f..%2.0f GeV", dR_->GetName(), decayMode_string.Data(), ptMin, ptMax);
      RooPlot* frame = dR_->frame(RooFit::Title(frameTitle.Data()), RooFit::Bins(50));
    
      datahist->plotOn(frame);

      prefitModel_[iPtBin]->plotOn(frame);
      prefitModel_[iPtBin]->plotOn(frame, RooFit::Components(*prefitLandau_[iPtBin]), RooFit::LineStyle(kDashed));
      prefitModel_[iPtBin]->plotOn(frame, RooFit::Components(*prefitGaussian_[iPtBin]), RooFit::LineStyle(kDashDotted));

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
      outputFileName_i = outputFileName_i.ReplaceAll(".", TString("_").Append(ptBinName).Append("."));
      
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
      Form("LandauScale_%s_%s_%s_%s", decayMode_string.Data(), "AllPt", dR_->GetName(), label_.Data());
    fitLandauScale_ = 
      new RooFormulaVar(LandauScaleName.Data(), LandauScaleName.Data(), "@0/@1", RooArgSet(*dR_, *fitParamLandauScaleFactor_));
    TString LandauName =
      Form("Landau_%s_%s_%s_%s", decayMode_string.Data(), "AllPt", dR_->GetName(), label_.Data());
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
      Form("Gaussian_%s_%s_%s_%s", decayMode_string.Data(), "AllPt", dR_->GetName(), label_.Data());
    fitGaussian_ =
      new RooGaussian(GaussianName.Data(), GaussianName.Data(), *dR_,
		      *fitParamGaussianMean_, *fitParamGaussianSigma_);
    
    fitParamMix_ = 
      buildMomDependentFitParameter("mix", fitParamCoeffMix_,
				    prefitCoeffValMix_, prefitCoeffErrMix_, true);
			
    TString modelName = Form("pdf_%s_%s_%s_%s", decayMode_string.Data(), "AllPt", dR_->GetName(), label_.Data());
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

    for ( int iPtBin = 0; iPtBin < (ptBinning_.GetSize() - 1); ++iPtBin ) {
      double ptMin = ptBinning_[iPtBin];
      double ptMax = ptBinning_[iPtBin + 1];

      if ( (iPtBin % 6) == 0 ) {
	canvas->Clear();
	canvas->Divide(2,3);
      }

      TVirtualPad* pad = canvas->cd((iPtBin % 6) + 1);
      pad->SetLogy();

      TString frameTitle = Form("%s %s: P_{T} = %2.0f..%2.0f GeV", dR_->GetName(), decayMode_string.Data(), ptMin, ptMax);
      RooPlot* frame = dR_->frame(RooFit::Title(frameTitle.Data()), RooFit::Bins(50));
      
      TString ptSliceName = Form("ptSlice%2.0fto%2.0f", ptMin, ptMax);
      std::stringstream ptSliceCut;
      ptSliceCut << ptMin << " < " << mom_->GetName() << " && " << mom_->GetName() << " < " << ptMax;
      
      RooAbsData* ptSliceDataset = dataset->reduce(ptSliceCut.str().data());
      std::cout << "Selected " << ptSliceDataset->numEntries() << " TTree entries" 
		<< " in slice " << Form("pT = %2.0f..%2.0f GeV", ptMin, ptMax) << "..." << std::endl;
      RooDataHist ptSliceDataset_binned("ptSliceDataset_binned", "ptSliceDataset_binned", RooArgSet(*mom_, *dR_), *ptSliceDataset);
      ptSliceDataset_binned.plotOn(frame);
     
      fitModel_->plotOn(frame, RooFit::ProjWData(RooArgSet(*mom_), ptSliceDataset_binned));
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

      if ( (iPtBin % 6) == 5 ) {
	canvas->Update();
	TString outputFileName_i = outputFileName;
	outputFileName_i = outputFileName_i.ReplaceAll(".", Form("_page%u.", (iPtBin / 6) + 1));
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
  TArrayD ptBinning_;

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

  //unsigned numPtBins = 12;
  //const double ptBinning[] = { 15., 20., 25., 30., 35., 40., 50., 60., 80., 100., 120., 160., 200. };
  unsigned numPtBins = 36;
  const double ptBinning[] = { 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75.,
			       80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135., 140., 
			       145., 150., 155., 160., 165., 170., 175., 180., 185., 190., 195. };

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

  RooRealVar leg1Energy("leg1Energy", "leg1Energy", 0., 300.);
  RooRealVar leg1Pt("leg1Pt", "leg1Pt", 0., 200.);
  RooRealVar leg1VisPt("leg1VisPt", "leg1VisPt", 0., 200.);
  RooRealVar leg1VisEta("leg1VisEta", "leg1VisEta", -2.5, +2.5);
  RooRealVar leg1DecayMode("leg1DecayMode", "leg1DecayMode", -1.5, +11.5);
  RooRealVar leg1VisInvisAngleLab("leg1VisInvisAngleLab", "leg1VisInvisAngleLab", 0., 0.50);
  RooRealVar leg1VisInvisDeltaRLab("leg1VisInvisDeltaRLab", "leg1VisInvisDeltaRLab", 0., 0.50);

  leg1Energy.setBins(100);
  leg1Pt.setBins(100);
  leg1VisInvisAngleLab.setBins(100);
  leg1VisInvisDeltaRLab.setBins(100);

  RooRealVar leg2Energy("leg2Energy", "leg2Energy", 0., 300.);
  RooRealVar leg2Pt("leg2Pt", "leg2Pt", 0., 200.);
  RooRealVar leg2VisPt("leg2VisPt", "leg2VisPt", 0., 200.);
  RooRealVar leg2VisEta("leg2VisEta", "leg2VisEta", -2.5, +2.5);
  RooRealVar leg2DecayMode("leg2DecayMode", "leg2DecayMode", -1.5, +11.5);
  RooRealVar leg2VisInvisAngleLab("leg2VisInvisAngleLab", "leg2VisInvisAngleLab", 0., 0.50);
  RooRealVar leg2VisInvisDeltaRLab("leg2VisInvisDeltaRLab", "leg2VisInvisDeltaRLab", 0., 0.50);

  leg2Energy.setBins(100);
  leg2Pt.setBins(100);
  leg2VisInvisAngleLab.setBins(100);
  leg2VisInvisDeltaRLab.setBins(100);

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

  RooAbsData* dataset_all = new RooDataSet("dataset", "datasetset", RooArgSet(variables), RooFit::Import(*dataTree));
  std::cout << "Processing " << dataset_all->numEntries() << " TTree entries..." << std::endl;
  
  RooAbsData* dataset_notNaN = dataset_all->reduce(nanFilter);
  std::cout << "--> Passing anti-NaN filter = " << dataset_notNaN->numEntries() << std::endl;
  
  RooAbsData* dataset_selected = 0;
  if ( runSelected ) {
    dataset_selected = dataset_notNaN->reduce(visMomCuts);
    std::cout << "--> Passing vis. Momentum Cuts = " << dataset_selected->numEntries() << std::endl;
  }

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
      RooAbsData* dataset_leg1_all = dataset_notNaN->reduce(leg1DecayModeSelection.Data());
      
      TString fitManagerName_leg1_dR_all = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg1", "dR", "all");
      if ( fitResults.find(fitManagerName_leg1_dR_all.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg1_dR_all.Data()] = 
	  new fitManager(&leg1Pt, &leg1VisInvisDeltaRLab, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg1", "dR", "all"), numPtBins, ptBinning);
      }
      fitManager* fitManager_leg1_dR_all = fitResults[fitManagerName_leg1_dR_all.Data()];
      TString outputFileName_leg1_dR_all = Form("plots/fitTauDecayKinePlots_%s_leg1_dR_all.eps", decayMode_string.Data());
      fitManager_leg1_dR_all->runPrefit(inputFileName_histograms, 
					"VisInvisDeltaRLab", "all", outputFileName_leg1_dR_all);
      fitManager_leg1_dR_all->runFit(dataset_leg1_all, outputFileName_leg1_dR_all);
      
      TString fitManagerName_leg1_angle_all = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg1", "angle", "all");
      if ( fitResults.find(fitManagerName_leg1_angle_all.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg1_angle_all.Data()] = 
	  new fitManager(&leg1Energy, &leg1VisInvisAngleLab, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg1", "angle", "all"), numPtBins, ptBinning);
      }
      fitManager* fitManager_leg1_angle_all = fitResults[fitManagerName_leg1_angle_all.Data()];
      TString outputFileName_leg1_angle_all = 
	Form("plots/fitTauDecayKinePlots_%s_leg1_angle_all.eps", decayMode_string.Data());
      fitManager_leg1_angle_all->runPrefit(inputFileName_histograms, 
					   "VisInvisAngleLab", "all", outputFileName_leg1_angle_all);
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
			 Form("%s_%s_%s", "leg1", "dR", "selected2"), numPtBins, ptBinning);
      }
      fitManager* fitManager_leg1_dR_selected = fitResults[fitManagerName_leg1_dR_selected.Data()];
      TString outputFileName_leg1_dR_selected = Form("plots/fitTauDecayKinePlots_%s_leg1_dR_selected.eps", decayMode_string.Data());
      fitManager_leg1_dR_selected->runPrefit(inputFileName_histograms, 
					   "VisInvisDeltaRLab", "selected2", outputFileName_leg1_dR_selected);
      fitManager_leg1_dR_selected->runFit(dataset_leg1_selected, outputFileName_leg1_dR_selected);
      
      TString fitManagerName_leg1_angle_selected = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg1", "angle", "selected");
      if ( fitResults.find(fitManagerName_leg1_angle_selected.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg1_angle_selected.Data()] = 
	  new fitManager(&leg1Energy, &leg1VisInvisAngleLab, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg1", "angle", "selected2"), numPtBins, ptBinning);
      }
      fitManager* fitManager_leg1_angle_selected = fitResults[fitManagerName_leg1_angle_selected.Data()];
      TString outputFileName_leg1_angle_selected = 
	Form("plots/fitTauDecayKinePlots_%s_leg1_angle_selected.eps", decayMode_string.Data());
      fitManager_leg1_angle_selected->runPrefit(inputFileName_histograms, 
						"VisInvisAngleLab", "selected2", outputFileName_leg1_angle_selected);
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
			 Form("%s_%s_%s", "leg2", "dR", "all"), numPtBins, ptBinning);
      }
      fitManager* fitManager_leg2_dR_all = fitResults[fitManagerName_leg2_dR_all.Data()];
      TString outputFileName_leg2_dR_all = Form("plots/fitTauDecayKinePlots_%s_leg2_dR_all.eps", decayMode_string.Data());
      fitManager_leg2_dR_all->runPrefit(inputFileName_histograms, 
					"VisInvisDeltaRLab", "all", outputFileName_leg2_dR_all);
      fitManager_leg2_dR_all->runFit(dataset_leg2_all, outputFileName_leg2_dR_all);
      
      TString fitManagerName_leg2_angle_all = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg2", "angle", "all");
      if ( fitResults.find(fitManagerName_leg2_angle_all.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg2_angle_all.Data()] = 
	  new fitManager(&leg2Energy, &leg2VisInvisAngleLab, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg2", "angle", "all"), numPtBins, ptBinning);
      }
      fitManager* fitManager_leg2_angle_all = fitResults[fitManagerName_leg2_angle_all.Data()];
      TString outputFileName_leg2_angle_all = 
	Form("plots/fitTauDecayKinePlots_%s_leg2_angle_all.eps", decayMode_string.Data());
      fitManager_leg2_angle_all->runPrefit(inputFileName_histograms, 
					   "VisInvisAngleLab", "all", outputFileName_leg2_angle_all);
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
			 Form("%s_%s_%s", "leg2", "dR", "selected2"), numPtBins, ptBinning);
      }
      fitManager* fitManager_leg2_dR_selected = fitResults[fitManagerName_leg2_dR_selected.Data()];
      TString outputFileName_leg2_dR_selected = Form("plots/fitTauDecayKinePlots_%s_leg2_dR_selected.eps", decayMode_string.Data());
      fitManager_leg2_dR_selected->runPrefit(inputFileName_histograms, 
					     "VisInvisDeltaRLab", "selected2", outputFileName_leg2_dR_selected);
      fitManager_leg2_dR_selected->runFit(dataset_leg2_selected, outputFileName_leg2_dR_selected);
      
      TString fitManagerName_leg2_angle_selected = Form("%s_%s_%s_%s", decayMode_string.Data(), "leg2", "angle", "selected");
      if ( fitResults.find(fitManagerName_leg2_angle_selected.Data()) == fitResults.end() ) {
	fitResults[fitManagerName_leg2_angle_selected.Data()] = 
	  new fitManager(&leg2Energy, &leg2VisInvisAngleLab, (*decayModeToRun), 
			 Form("%s_%s_%s", "leg2", "angle", "selected2"), numPtBins, ptBinning);
      }
      fitManager* fitManager_leg2_angle_selected = fitResults[fitManagerName_leg2_angle_selected.Data()];
      TString outputFileName_leg2_angle_selected = 
	Form("plots/fitTauDecayKinePlots_%s_leg2_angle_selected.eps", decayMode_string.Data());
      fitManager_leg2_angle_selected->runPrefit(inputFileName_histograms, 
						"VisInvisAngleLab", "selected2", outputFileName_leg2_angle_selected);
      fitManager_leg2_angle_selected->runFit(dataset_leg2_selected, outputFileName_leg2_angle_selected);
    }

//--- write results of prefit to ROOT file
    TString outputFileName_prefit = Form("fitTauDecayKinePlots_%s%s.root", decayMode_string.Data(), selection_string.Data());
    TFile* outputFile_prefit = new TFile(outputFileName_prefit.Data(), "RECREATE");
    TDirectory* dirVisInvisDeltaRLab = outputFile_prefit->mkdir("VisInvisDeltaRLab");  
    TDirectory* dirVisInvisDeltaRLab_all = dirVisInvisDeltaRLab->mkdir("all");
    TDirectory* dirVisInvisDeltaRLab_selected = dirVisInvisDeltaRLab->mkdir("selected");
    TDirectory* dirVisInvisAngleLab = outputFile_prefit->mkdir("VisInvisAngleLab");
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
