
#include "TauAnalysis/FittingTools/interface/SmoothLandau_x_GaussPdf.h"

#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooGenericPdf.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooCmdArg.h"
#include "RooFit.h"
#include "RooLinkedList.h"
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

#include <iostream>
#include <iomanip>
#include <string>

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

//
//-------------------------------------------------------------------------------
//

void storePrefitResults(RooRealVar* fitParameter, double ptMin, double ptMax, TGraphErrors* graph, double* mixFactor = 0)
{
  std::cout << "<storePrefitResults>:" << std::endl;

  int iPoint = graph->GetN();

  double pt = 0.5*(ptMin + ptMax);

  graph->SetPoint(iPoint, pt, fitParameter->getVal());

//--- CV: enlarge uncertainties by factor 10 and weight coefficients of smearedLandau/skewedGaussian 
//        by corresponding "mix" fraction, in order to improve convergence when fitting TGraph by TF1 
  double err = fitParameter->getError();
  err *= 10.;
  if ( mixFactor ) err *= (*mixFactor);

  graph->SetPointError(iPoint, 0.5*TMath::Abs(ptMax - ptMin), err);
}

void storeFitResults(std::vector<RooRealVar*>& fitParamCoeff, std::vector<double>& fitCoeffVal, std::vector<double>& fitCoeffErr)
{
  std::cout << "<storeFitResults>:" << std::endl;

  unsigned numFitParameter = fitParamCoeff.size();
  std::cout << " numFitParameter = " << numFitParameter << std::endl;
  
  fitCoeffVal.resize(numFitParameter);
  fitCoeffErr.resize(numFitParameter);
 
  for ( unsigned iParameter = 0; iParameter < numFitParameter; ++iParameter ) {
    std::cout << "iParameter = " << iParameter << ": fitParamCoeff[iParameter] = " << fitParamCoeff[iParameter] << std::endl;
    fitCoeffVal[iParameter] = fitParamCoeff[iParameter]->getVal();
    fitCoeffErr[iParameter] = fitParamCoeff[iParameter]->getError();
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
      prefitParamSmearedLandauMP_(numPtBins),
      prefitParamSmearedLandauWidth_(numPtBins),
      prefitSmearedLandauPdf_(numPtBins), 
      prefitParamSmearedLandauGMean_(numPtBins), 
      prefitParamSmearedLandauGSigma_(numPtBins), 
      prefitSmearedLandauGPdf_(numPtBins), 
      prefitSmearedLandau_(numPtBins), 
      prefitParamSkewedGaussianMean_(numPtBins),
      prefitParamSkewedGaussianSigma_(numPtBins),
      prefitParamSkewedGaussianAlpha_(numPtBins),
      prefitSkewedGaussian_(numPtBins), 
      prefitParamMix_(numPtBins),
      prefitModel_(numPtBins),
      fitParamCoeffSmearedLandauMP_(6),
      fitParamSmearedLandauMP_(0),
      fitParamCoeffSmearedLandauWidth_(6),
      fitParamSmearedLandauWidth_(0),     
      fitSmearedLandauPdf_(0),
      fitParamSmearedLandauGMean_(0),
      fitParamCoeffSmearedLandauGSigma_(6),
      fitParamSmearedLandauGSigma_(0),
      fitSmearedLandauGPdf_(0),
      fitSmearedLandau_(0),
      fitParamCoeffSkewedGaussianMean_(6),
      fitParamSkewedGaussianMean_(0),
      fitParamCoeffSkewedGaussianSigma_(6),
      fitParamSkewedGaussianSigma_(0),
      fitParamCoeffSkewedGaussianAlpha_(6),
      fitParamSkewedGaussianAlpha_(0),
      fitSkewedGaussian_(0),
      fitParamCoeffMix_(6),
      fitParamMix_(0),
      fitModel_(0)
  {
    std::cout << "<fitManager::fitManager>:" << std::endl;
    std::cout << " mom = " << mom_->GetName() << std::endl;
    std::cout << " dR = " << dR_->GetName() << std::endl;
    std::cout << " decayMode = " << getDecayMode_string(decayMode_) << std::endl;

    prefitResultSmearedLandauMP_ = new TGraphErrors();
    prefitResultSmearedLandauWidth_ = new TGraphErrors();
    prefitResultSmearedLandauGSigma_ = new TGraphErrors();
    prefitResultSkewedGaussianMean_ = new TGraphErrors();
    prefitResultSkewedGaussianSigma_ = new TGraphErrors();
    prefitResultSkewedGaussianAlpha_ = new TGraphErrors();
    prefitResultMix_ = new TGraphErrors();

//--- parametrize pt/energy dependence of fitParameters by (orthogonal) Chebyshev polynomials
//   ( cf. http://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html )
    TString momParametrizationTFormula = 
      "(1./(x*x*x*x))*([0] + [1]*x + [2]*(2.0*x*x - 1.0)" 
      " + [3]*(4.0*x*x*x - 3.0*x) + [4]*(8.0*x*x*x*x - 8.0*x*x + 1.0))";
    momParametrization_forTFormula_ = momParametrizationTFormula;
    momParametrizationMix_forTFormula_ = TString("0.5*(1.0 + TMath::TanH(").Append(momParametrization_forTFormula_).Append("))");
    TString momParametrizationRooFit = 
      "(1./(@5*@5*@5*@5))*(@0 + @1*@5 + @2*(2.0*@5*@5 - 1.0)" 
      " + @3*(4.0*@5*@5*@5 - 3.0*@5) + @4*(8.0*@5*@5*@5*@5 - 8.0*@5*@5 + 1.0))";
    momParametrization_forRooFit_ = momParametrizationRooFit;
    momParametrizationMix_forRooFit_ = TString("0.5*(1.0 + TMath::TanH(").Append(momParametrization_forRooFit_).Append("))");
    
    momParamSmearedLandauMP_ = bookTF1("smearedLandauMP", momParametrization_forTFormula_);
    momParamSmearedLandauWidth_ = bookTF1("smearedLandauWidth", momParametrization_forTFormula_);
    momParamSmearedLandauGSigma_ = bookTF1("smearedLandauGSigma", momParametrization_forTFormula_);
    momParamSkewedGaussianMean_ = bookTF1("skewedGaussian_mean", momParametrization_forTFormula_);
    momParamSkewedGaussianSigma_ = bookTF1("skewedGaussian_sigma", momParametrization_forTFormula_);
    momParamSkewedGaussianAlpha_ = bookTF1("skewedGaussian_alpha", momParametrization_forTFormula_);
    momParamMix_ = bookTF1("mix", momParametrizationMix_forTFormula_);
  }

  ~fitManager()
  {
    clearCollection(prefitParamSmearedLandauMP_);
    clearCollection(prefitParamSmearedLandauWidth_);
    clearCollection(prefitSmearedLandauPdf_);
    clearCollection(prefitParamSmearedLandauGMean_);
    clearCollection(prefitParamSmearedLandauGSigma_);
    clearCollection(prefitSmearedLandauGPdf_);
    clearCollection(prefitSmearedLandau_);
    clearCollection(prefitParamSkewedGaussianMean_);
    clearCollection(prefitParamSkewedGaussianSigma_);
    clearCollection(prefitParamSkewedGaussianAlpha_);
    clearCollection(prefitSkewedGaussian_);
    clearCollection(prefitParamMix_);
    clearCollection(prefitModel_);

    delete prefitResultSmearedLandauMP_;
    delete prefitResultSmearedLandauWidth_;
    delete prefitResultSmearedLandauGSigma_;
    delete prefitResultSkewedGaussianMean_;
    delete prefitResultSkewedGaussianSigma_;
    delete prefitResultSkewedGaussianAlpha_;

    delete momParamSmearedLandauMP_;
    delete momParamSmearedLandauWidth_;
    delete momParamSmearedLandauGSigma_;
    delete momParamSkewedGaussianMean_;
    delete momParamSkewedGaussianSigma_;
    delete momParamSkewedGaussianAlpha_;
    delete momParamMix_;
    
    clearCollection(fitParamCoeffSmearedLandauMP_);
    delete fitParamSmearedLandauMP_;
    clearCollection(fitParamCoeffSmearedLandauWidth_);
    delete fitParamSmearedLandauWidth_;
    delete fitSmearedLandauPdf_;
    delete fitParamSmearedLandauGMean_;
    clearCollection(fitParamCoeffSmearedLandauGSigma_);
    delete fitParamSmearedLandauGSigma_;
    delete fitSmearedLandauGPdf_;
    delete fitSmearedLandau_;
    clearCollection(fitParamCoeffSkewedGaussianMean_);
    delete fitParamSkewedGaussianMean_;
    clearCollection(fitParamCoeffSkewedGaussianSigma_);
    delete fitParamSkewedGaussianSigma_;
    clearCollection(fitParamCoeffSkewedGaussianAlpha_);
    delete fitParamSkewedGaussianAlpha_;
    delete fitSkewedGaussian_;
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

      prefitParamSmearedLandauMP_[iPtBin] = buildPrefitParameter("smearedLandauMP", ptBinName, 1.2*histogramMean, -0.20, +0.20);
      prefitParamSmearedLandauWidth_[iPtBin] = buildPrefitParameter("smearedLandauWidth", ptBinName, histogramRMS, 0., 0.20);
      TString smearedLandauPdfName =
	Form("smearedLandauPdf_%s_%s_%s_%s", decayMode_string.Data(), ptBinName.Data(), dR_->GetName(), label_.Data());
      prefitSmearedLandauPdf_[iPtBin] = 
	new RooLandau(smearedLandauPdfName.Data(), smearedLandauPdfName.Data(), *dR_, 
		      *prefitParamSmearedLandauMP_[iPtBin], *prefitParamSmearedLandauWidth_[iPtBin]);
      prefitParamSmearedLandauGMean_[iPtBin] = buildConstPrefitParameter("smearedLandauGMean", ptBinName, 0.);
      prefitParamSmearedLandauGSigma_[iPtBin] = buildPrefitParameter("smearedLandauGSigma", ptBinName, histogramRMS, 0., 0.20);
      TString smearedLandauGPdfName =
	Form("smearedLandauGPdf_%s_%s_%s_%s", decayMode_string.Data(), ptBinName.Data(), dR_->GetName(), label_.Data());
      prefitSmearedLandauGPdf_[iPtBin] = 
	new RooGaussian(smearedLandauGPdfName.Data(), smearedLandauGPdfName.Data(), *dR_, 
			*prefitParamSmearedLandauGMean_[iPtBin], *prefitParamSmearedLandauGSigma_[iPtBin]);
      TString smearedLandauName =
	Form("smearedLandau_%s_%s_%s_%s", decayMode_string.Data(), ptBinName.Data(), dR_->GetName(), label_.Data());
      prefitSmearedLandau_[iPtBin] = 
	new RooFFTConvPdf(smearedLandauName.Data(), smearedLandauName.Data(), *dR_, 
			  *prefitSmearedLandauPdf_[iPtBin], *prefitSmearedLandauGPdf_[iPtBin]);

      prefitParamSkewedGaussianMean_[iPtBin] = buildPrefitParameter("skewedGaussianMean", ptBinName, 0.8*histogramMean, 0., 0.20);
      prefitParamSkewedGaussianSigma_[iPtBin] = buildPrefitParameter("skewedGaussianSigma", ptBinName, histogramRMS, 0.01, 0.20);
      prefitParamSkewedGaussianAlpha_[iPtBin] = buildPrefitParameter("skewedGaussianAlpha", ptBinName, 0., -2., +2.);
      TString skewedGaussianName = 
	Form("skewedGaussian_%s_%s_%s_%s", decayMode_string.Data(), ptBinName.Data(), dR_->GetName(), label_.Data());
      prefitSkewedGaussian_[iPtBin] = 
	new RooGenericPdf(skewedGaussianName.Data(), skewedGaussianName.Data(), 
			  "TMath::Gaus(@0, @1, @2)*(1.0 + TMath::Erf(@3*((@0 - @1)/@2)/TMath::Sqrt(2)))",
			  RooArgList(*dR_, *prefitParamSkewedGaussianMean_[iPtBin], *prefitParamSkewedGaussianSigma_[iPtBin], 
				     *prefitParamSkewedGaussianAlpha_[iPtBin])); 
      
      prefitParamMix_[iPtBin] = buildPrefitParameter("mix", ptBinName, 0.50, 0., 1.);
  
      TString modelName = Form("pdf_%s_%s_%s_%s", decayMode_string.Data(), ptBinName.Data(), dR_->GetName(), label_.Data());
      prefitModel_[iPtBin] = 
	new RooAddPdf(modelName.Data(), modelName.Data(), 
		      *prefitSmearedLandau_[iPtBin], *prefitSkewedGaussian_[iPtBin], *prefitParamMix_[iPtBin]);
  
      RooConstVar* skewedGaussianAlphaConstraint_value =
        new RooConstVar("skewedGaussianAlphaConstraint_value", "skewedGaussianAlphaConstraint_value", 0.);
      RooConstVar* skewedGaussianAlphaConstraint_sigma =
        new RooConstVar("skewedGaussianAlphaConstraint_sigma", "skewedGaussianAlphaConstraint_sigma", 1.);
      RooGaussian* skewedGaussianAlphaConstraint_pdf =
        new RooGaussian("skewedGaussianAlphaConstraint_pdf", "skewedGaussianAlphaConstraint_pdf",
			*prefitParamSkewedGaussianAlpha_[iPtBin], 
			*skewedGaussianAlphaConstraint_value, *skewedGaussianAlphaConstraint_sigma);
    
      RooLinkedList options;
      options.Add(new RooCmdArg(RooFit::PrintLevel(-1)));
      options.Add(new RooCmdArg(RooFit::PrintEvalErrors(-1)));
      options.Add(new RooCmdArg(RooFit::Warnings(-1)));
      options.Add(new RooCmdArg(RooFit::ExternalConstraints(RooArgSet(*skewedGaussianAlphaConstraint_pdf))));
      
      prefitModel_[iPtBin]->fitTo(*datahist, options); 

      printPrefitParameter("smearedLandauMP", prefitParamSmearedLandauMP_[iPtBin]);
      printPrefitParameter("smearedLandauWidth", prefitParamSmearedLandauWidth_[iPtBin]);
      printPrefitParameter("smearedLandauGSigma", prefitParamSmearedLandauGSigma_[iPtBin]);
      printPrefitParameter("skewedGaussianMean", prefitParamSkewedGaussianMean_[iPtBin]);
      printPrefitParameter("skewedGaussianSigma", prefitParamSkewedGaussianSigma_[iPtBin]);
      printPrefitParameter("skewedGaussianAlpha", prefitParamSkewedGaussianAlpha_[iPtBin]);
      printPrefitParameter("mix", prefitParamMix_[iPtBin]);

      double mixSmearedLandau = prefitParamMix_[iPtBin]->getVal();
      double scaleSmearedLandau = ( mixSmearedLandau > 0. ) ? 1./mixSmearedLandau : 100.;
      storePrefitResults(prefitParamSmearedLandauMP_[iPtBin], ptMin, ptMax, prefitResultSmearedLandauMP_, &scaleSmearedLandau);
      storePrefitResults(prefitParamSmearedLandauWidth_[iPtBin], ptMin, ptMax, prefitResultSmearedLandauWidth_, &scaleSmearedLandau);
      storePrefitResults(prefitParamSmearedLandauGSigma_[iPtBin], ptMin, ptMax, prefitResultSmearedLandauGSigma_, &scaleSmearedLandau);
      double mixSkewedGaussian = 1. - mixSmearedLandau;
      double scaleSkewedGaussian = ( mixSkewedGaussian > 0. ) ? 1./mixSkewedGaussian : 100.;
      storePrefitResults(prefitParamSkewedGaussianMean_[iPtBin], ptMin, ptMax, prefitResultSkewedGaussianMean_, &scaleSkewedGaussian);
      storePrefitResults(prefitParamSkewedGaussianSigma_[iPtBin], ptMin, ptMax, prefitResultSkewedGaussianSigma_, &scaleSkewedGaussian);
      storePrefitResults(prefitParamSkewedGaussianAlpha_[iPtBin], ptMin, ptMax, prefitResultSkewedGaussianAlpha_, &scaleSkewedGaussian);
      storePrefitResults(prefitParamMix_[iPtBin], ptMin, ptMax, prefitResultMix_);
         
      canvas->Clear();
      canvas->SetLogy();
    
      TString frameTitle = Form("%s %s: P_{T} = %2.0f..%2.0f GeV", dR_->GetName(), decayMode_string.Data(), ptMin, ptMax);
      RooPlot* frame = dR_->frame(RooFit::Title(frameTitle.Data()), RooFit::Bins(50));
    
      datahist->plotOn(frame);

      prefitModel_[iPtBin]->plotOn(frame);
      prefitModel_[iPtBin]->plotOn(frame, RooFit::Components(*prefitSmearedLandau_[iPtBin]), RooFit::LineStyle(kDashed));
      prefitModel_[iPtBin]->plotOn(frame, RooFit::Components(*prefitSkewedGaussian_[iPtBin]), RooFit::LineStyle(kDashDotted));

      RooConstVar* approxSmearedLandauArea = buildConstPrefitParameter("approxSmearedLandauArea", ptBinName, 1.);
      SmoothLandau_x_GaussPdf* approxSmearedLandau = 
	new SmoothLandau_x_GaussPdf("approxSmearedLandau", "approxSmearedLandau", *dR_,
				    *prefitParamSmearedLandauWidth_[iPtBin], *prefitParamSmearedLandauMP_[iPtBin], 
				    *approxSmearedLandauArea, *prefitParamSmearedLandauGSigma_[iPtBin]);
      approxSmearedLandau->plotOn(frame, RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));

      frame->Draw();

      TLegend legend(0.60, 0.64, 0.89, 0.89, "", "brNDC");
      legend.SetBorderSize(0);
      legend.SetFillColor(0);
      TObject* graphSmearedLandau        = 0;
      TObject* graphSmearedLandau_approx = 0;
      TObject* graphSkewedGaussian       = 0;
      TObject* graphSum                  = 0;
      TObject* graphData                 = 0;
      for ( int iItem = 0; iItem < frame->numItems(); ++iItem ) {
	TString itemName = frame->nameOf(iItem);
	TObject* item = frame->findObject(itemName.Data());
	if      ( itemName.Contains("smearedLandau")       ) graphSmearedLandau        = item;
	else if ( itemName.Contains("approxSmearedLandau") ) graphSmearedLandau_approx = item;
	else if ( itemName.Contains("skewedGaussian")      ) graphSkewedGaussian       = item;
	else if ( itemName.Contains("pdf")                 ) graphSum                  = item;
	else if ( itemName.Contains("datahist")            ) graphData                 = item; 
      }
      if ( graphSmearedLandau        != 0 ) legend.AddEntry(graphSmearedLandau,        "Landau #otimes Gaussian",         "l");
      if ( graphSmearedLandau_approx != 0 ) legend.AddEntry(graphSmearedLandau_approx, "Landau #otimes Gaussian approx.", "l");
      if ( graphSkewedGaussian       != 0 ) legend.AddEntry(graphSkewedGaussian,       "skewed Gaussian",                 "l");
      if ( graphSum                  != 0 ) legend.AddEntry(graphSum,                  "Sum",                             "l");
      if ( graphData                 != 0 ) legend.AddEntry(graphData,                 "TAUOLA",                          "p");
      legend.Draw();

      TString outputFileName_i = outputFileName;
      outputFileName_i = outputFileName_i.ReplaceAll(".", TString("_").Append(ptBinName).Append("."));
      
      canvas->Update();
  
      canvas->SaveAs(outputFileName_i.Data());

      delete approxSmearedLandauArea;
      delete approxSmearedLandau;

      delete datahist;
    }

    delete canvas;

    delete inputFile;
    
    fitTF1(momParamSmearedLandauMP_, prefitResultSmearedLandauMP_, 
	   prefitCoeffValSmearedLandauMP_, prefitCoeffErrSmearedLandauMP_);
    fitTF1(momParamSmearedLandauWidth_, prefitResultSmearedLandauWidth_, 
	   prefitCoeffValSmearedLandauWidth_, prefitCoeffErrSmearedLandauWidth_);
    fitTF1(momParamSmearedLandauGSigma_, prefitResultSmearedLandauGSigma_, 
	   prefitCoeffValSmearedLandauGSigma_, prefitCoeffErrSmearedLandauGSigma_);
    fitTF1(momParamSkewedGaussianMean_, prefitResultSkewedGaussianMean_, 
	   prefitCoeffValSkewedGaussianMean_, prefitCoeffErrSkewedGaussianMean_);
    fitTF1(momParamSkewedGaussianSigma_, prefitResultSkewedGaussianSigma_, 
	   prefitCoeffValSkewedGaussianSigma_, prefitCoeffErrSkewedGaussianSigma_);
    fitTF1(momParamSkewedGaussianAlpha_, prefitResultSkewedGaussianAlpha_, 
	   prefitCoeffValSkewedGaussianAlpha_, prefitCoeffErrSkewedGaussianAlpha_);
    fitTF1(momParamMix_, prefitResultMix_, 
	   prefitCoeffValMix_, prefitCoeffErrMix_);
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
		<< " startValue = " << prefitParamCoeffValues[iCoeff] << std::endl;
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

    fitParamSmearedLandauMP_ = 
      buildMomDependentFitParameter("smearedLandauMP", fitParamCoeffSmearedLandauMP_,
				    prefitCoeffValSmearedLandauMP_, prefitCoeffErrSmearedLandauMP_);
    fitParamSmearedLandauWidth_ = 
      buildMomDependentFitParameter("smearedLandauWidth", fitParamCoeffSmearedLandauWidth_,
				    prefitCoeffValSmearedLandauWidth_, prefitCoeffErrSmearedLandauWidth_);
    TString smearedLandauPdfName =
      Form("smearedLandauPdf_%s_%s_%s_%s", decayMode_string.Data(), "AllPt", dR_->GetName(), label_.Data());
    fitSmearedLandauPdf_ = 
      new RooLandau(smearedLandauPdfName.Data(), smearedLandauPdfName.Data(), *dR_, 
		    *fitParamSmearedLandauMP_, *fitParamSmearedLandauWidth_);
    fitParamSmearedLandauGMean_ = buildConstPrefitParameter("smearedLandauGMean", "AllPt", 0.);
    fitParamSmearedLandauGSigma_ = 
      buildMomDependentFitParameter("smearedLandauGSigma", fitParamCoeffSmearedLandauGSigma_,
				    prefitCoeffValSmearedLandauGSigma_, prefitCoeffErrSmearedLandauGSigma_);
    TString smearedLandauGPdfName =
      Form("smearedLandauGPdf_%s_%s_%s_%s", decayMode_string.Data(), "AllPt", dR_->GetName(), label_.Data());
    fitSmearedLandauGPdf_ = 
      new RooGaussian(smearedLandauGPdfName.Data(), smearedLandauGPdfName.Data(), *dR_, 
		      *fitParamSmearedLandauGMean_, *fitParamSmearedLandauGSigma_);
    TString smearedLandauName =
      Form("smearedLandau_%s_%s_%s_%s", decayMode_string.Data(), "AllPt", dR_->GetName(), label_.Data());
    fitSmearedLandau_ = 
      new RooFFTConvPdf(smearedLandauName.Data(), smearedLandauName.Data(), *dR_, 
			*fitSmearedLandauPdf_, *fitSmearedLandauGPdf_);
    
    fitParamSkewedGaussianMean_ = 
      buildMomDependentFitParameter("skewedGaussianMean", fitParamCoeffSkewedGaussianMean_,
				    prefitCoeffValSkewedGaussianMean_, prefitCoeffErrSkewedGaussianMean_);
    fitParamSkewedGaussianSigma_ = 
      buildMomDependentFitParameter("skewedGaussianSigma", fitParamCoeffSkewedGaussianSigma_,
				    prefitCoeffValSkewedGaussianSigma_, prefitCoeffErrSkewedGaussianSigma_);
    fitParamSkewedGaussianAlpha_ = 
      buildMomDependentFitParameter("skewedGaussianAlpha", fitParamCoeffSkewedGaussianAlpha_,
				    prefitCoeffValSkewedGaussianAlpha_, prefitCoeffErrSkewedGaussianAlpha_);
    TString skewedGaussianName = 
      Form("skewedGaussian_%s_%s_%s_%s", decayMode_string.Data(), "AllPt", dR_->GetName(), label_.Data());
    fitSkewedGaussian_ =
      new RooGenericPdf(skewedGaussianName.Data(), skewedGaussianName.Data(), 
			"TMath::Gaus(@0, @1, @2)*(1.0 + TMath::Erf(@3*((@0 - @1)/@2)/TMath::Sqrt(2)))",
			RooArgList(*dR_, *fitParamSkewedGaussianMean_, *fitParamSkewedGaussianSigma_, 
				   *fitParamSkewedGaussianAlpha_));
    
    fitParamMix_ = 
      buildMomDependentFitParameter("mix", fitParamCoeffMix_,
				    prefitCoeffValMix_, prefitCoeffErrMix_, true);
			
    TString modelName = Form("pdf_%s_%s_%s_%s", decayMode_string.Data(), "AllPt", dR_->GetName(), label_.Data());
    fitModel_ = 
      new RooAddPdf(modelName.Data(), modelName.Data(), 
		    *fitSmearedLandau_, *fitSkewedGaussian_, *fitParamMix_);
  
    RooLinkedList options;
    options.Add(new RooCmdArg(RooFit::ConditionalObservables(*mom_)));
    options.Add(new RooCmdArg(RooFit::PrintLevel(-1)));
    options.Add(new RooCmdArg(RooFit::PrintEvalErrors(-1)));
    options.Add(new RooCmdArg(RooFit::Warnings(-1)));
    
    fitModel_->fitTo(*dataset, options); 

    storeFitResults(fitParamCoeffSmearedLandauMP_, fitCoeffValSmearedLandauMP_, fitCoeffErrSmearedLandauMP_);
    storeFitResults(fitParamCoeffSmearedLandauWidth_, fitCoeffValSmearedLandauWidth_, fitCoeffErrSmearedLandauWidth_);
    storeFitResults(fitParamCoeffSmearedLandauGSigma_, fitCoeffValSmearedLandauGSigma_, fitCoeffErrSmearedLandauGSigma_);
    storeFitResults(fitParamCoeffSkewedGaussianMean_, fitCoeffValSkewedGaussianMean_, fitCoeffErrSkewedGaussianMean_);
    storeFitResults(fitParamCoeffSkewedGaussianSigma_, fitCoeffValSkewedGaussianSigma_, fitCoeffErrSkewedGaussianSigma_);
    storeFitResults(fitParamCoeffSkewedGaussianAlpha_, fitCoeffValSkewedGaussianAlpha_, fitCoeffErrSkewedGaussianAlpha_);
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
      
      TString ptBinName = Form("ptBin%2.0fto%2.0f", ptMin, ptMax);
      mom_->setRange(ptBinName.Data(), ptMin, ptMax);

      dataset->plotOn(frame, RooFit::CutRange(ptBinName));
      
      fitModel_->plotOn(frame, RooFit::ProjWData(dataset_binned), RooFit::ProjectionRange(ptBinName));
      
      frame->Draw();

      if ( (iPtBin % 6) == 5 ) {
	canvas->Update();
	TString outputFileName_i = outputFileName;
	outputFileName_i = outputFileName_i.ReplaceAll(".", Form("_page%u.", (iPtBin / 6) + 1));
	canvas->SaveAs(outputFileName_i.Data());
      }
    }

    delete canvas;
  }

  void writeFitResults(std::ostream& stream, const TString& psetName)
  {
    std::cout << "<fitManager::writeFitResults>:" << std::endl;

    stream << psetName << " = cms.PSet(";
    writeFitParameter(stream, "smearedLandauMP", momParametrization_forTFormula_, fitCoeffValSmearedLandauMP_);
    writeFitParameter(stream, "smearedLandauWidth", momParametrization_forTFormula_, fitCoeffValSmearedLandauWidth_);
    writeFitParameter(stream, "smearedLandauGSigma", momParametrization_forTFormula_, fitCoeffValSmearedLandauGSigma_);
    writeFitParameter(stream, "skewedGaussianMean", momParametrization_forTFormula_, fitCoeffValSkewedGaussianMean_);
    writeFitParameter(stream, "skewedGaussianSigma", momParametrization_forTFormula_, fitCoeffValSkewedGaussianSigma_);
    writeFitParameter(stream, "skewedGaussianAlpha", momParametrization_forTFormula_, fitCoeffErrSkewedGaussianAlpha_);
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

  std::vector<RooRealVar*> prefitParamSmearedLandauMP_;
  std::vector<RooRealVar*> prefitParamSmearedLandauWidth_;
  std::vector<RooLandau*> prefitSmearedLandauPdf_;
  std::vector<RooConstVar*> prefitParamSmearedLandauGMean_;
  std::vector<RooRealVar*> prefitParamSmearedLandauGSigma_;
  std::vector<RooGaussian*> prefitSmearedLandauGPdf_;
  std::vector<RooFFTConvPdf*> prefitSmearedLandau_;
  std::vector<RooRealVar*> prefitParamSkewedGaussianMean_;
  std::vector<RooRealVar*> prefitParamSkewedGaussianSigma_;
  std::vector<RooRealVar*> prefitParamSkewedGaussianAlpha_;
  std::vector<RooGenericPdf*> prefitSkewedGaussian_;
  std::vector<RooRealVar*> prefitParamMix_;
  std::vector<RooAddPdf*> prefitModel_;

  TGraphErrors* prefitResultSmearedLandauMP_;
  TGraphErrors* prefitResultSmearedLandauWidth_;
  TGraphErrors* prefitResultSmearedLandauGSigma_;
  TGraphErrors* prefitResultSkewedGaussianMean_;
  TGraphErrors* prefitResultSkewedGaussianSigma_;
  TGraphErrors* prefitResultSkewedGaussianAlpha_;
  TGraphErrors* prefitResultMix_;

  TString momParametrization_forTFormula_;
  TString momParametrizationMix_forTFormula_;
  TString momParametrization_forRooFit_;
  TString momParametrizationMix_forRooFit_;

  TF1* momParamSmearedLandauMP_;
  TF1* momParamSmearedLandauWidth_;
  TF1* momParamSmearedLandauGSigma_;
  TF1* momParamSkewedGaussianMean_;
  TF1* momParamSkewedGaussianSigma_;
  TF1* momParamSkewedGaussianAlpha_;
  TF1* momParamMix_;

  std::vector<double> prefitCoeffValSmearedLandauMP_;
  std::vector<double> prefitCoeffErrSmearedLandauMP_;
  std::vector<double> prefitCoeffValSmearedLandauWidth_;
  std::vector<double> prefitCoeffErrSmearedLandauWidth_;
  std::vector<double> prefitCoeffValSmearedLandauGSigma_;
  std::vector<double> prefitCoeffErrSmearedLandauGSigma_;
  std::vector<double> prefitCoeffValSkewedGaussianMean_;
  std::vector<double> prefitCoeffErrSkewedGaussianMean_;
  std::vector<double> prefitCoeffValSkewedGaussianSigma_;
  std::vector<double> prefitCoeffErrSkewedGaussianSigma_;
  std::vector<double> prefitCoeffValSkewedGaussianAlpha_;
  std::vector<double> prefitCoeffErrSkewedGaussianAlpha_;
  std::vector<double> prefitCoeffValMix_;
  std::vector<double> prefitCoeffErrMix_;

  
  std::vector<RooRealVar*> fitParamCoeffSmearedLandauMP_;
  RooAbsReal* fitParamSmearedLandauMP_;
  std::vector<RooRealVar*> fitParamCoeffSmearedLandauWidth_;
  RooAbsReal* fitParamSmearedLandauWidth_;
  RooLandau* fitSmearedLandauPdf_;
  RooConstVar* fitParamSmearedLandauGMean_;
  std::vector<RooRealVar*> fitParamCoeffSmearedLandauGSigma_;
  RooAbsReal* fitParamSmearedLandauGSigma_;
  RooGaussian* fitSmearedLandauGPdf_;
  RooFFTConvPdf* fitSmearedLandau_;
  std::vector<RooRealVar*> fitParamCoeffSkewedGaussianMean_;
  RooAbsReal* fitParamSkewedGaussianMean_;
  std::vector<RooRealVar*> fitParamCoeffSkewedGaussianSigma_;
  RooAbsReal* fitParamSkewedGaussianSigma_;
  std::vector<RooRealVar*> fitParamCoeffSkewedGaussianAlpha_;
  RooAbsReal* fitParamSkewedGaussianAlpha_;
  RooGenericPdf* fitSkewedGaussian_;
  std::vector<RooRealVar*> fitParamCoeffMix_;
  RooAbsReal*  fitParamMix_;
  RooAddPdf* fitModel_;

  std::vector<double> fitCoeffValSmearedLandauMP_;
  std::vector<double> fitCoeffErrSmearedLandauMP_;
  std::vector<double> fitCoeffValSmearedLandauWidth_;
  std::vector<double> fitCoeffErrSmearedLandauWidth_;
  std::vector<double> fitCoeffValSmearedLandauGSigma_;
  std::vector<double> fitCoeffErrSmearedLandauGSigma_;
  std::vector<double> fitCoeffValSkewedGaussianMean_;
  std::vector<double> fitCoeffErrSkewedGaussianMean_;
  std::vector<double> fitCoeffValSkewedGaussianSigma_;
  std::vector<double> fitCoeffErrSkewedGaussianSigma_;
  std::vector<double> fitCoeffValSkewedGaussianAlpha_;
  std::vector<double> fitCoeffErrSkewedGaussianAlpha_;
  std::vector<double> fitCoeffValMix_;
  std::vector<double> fitCoeffErrMix_;
};

//
//-------------------------------------------------------------------------------
//

int main(int argc, const char* argv[])
{
  //TString inputFileNames_ntuple = "/data2/friis/PtBalanceNtupleData_v2/ptBalanceData_*.root";
  TString inputFileNames_ntuple = "/data2/friis/PtBalanceNtupleData_v2/ptBalanceData_mass_200_ggAH.root";

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
  RooRealVar leg1VisInvisAngleLab("leg1VisInvisAngleLab", "leg1VisInvisAngleLab", 0., 0.50);
  RooRealVar leg1VisInvisDeltaRLab("leg1VisInvisDeltaRLab", "leg1VisInvisDeltaRLab", 0., 0.50);

  RooRealVar leg2Energy("leg2Energy", "leg2Energy", 0., 350.);
  RooRealVar leg2Pt("leg2Pt", "leg2Pt", 0., 250.);
  RooRealVar leg2VisPt("leg2VisPt", "leg2VisPt", 0., 200.);
  RooRealVar leg2VisEta("leg2VisEta", "leg2VisEta", -2.5, +2.5);
  RooRealVar leg2DecayMode("leg2DecayMode", "leg2DecayMode", -1.5, +11.5);
  RooRealVar leg2VisInvisAngleLab("leg2VisInvisAngleLab", "leg2VisInvisAngleLab", 0., 0.50);
  RooRealVar leg2VisInvisDeltaRLab("leg2VisInvisDeltaRLab", "leg2VisInvisDeltaRLab", 0., 0.50);

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

  //RooAbsData* dataset_selected = dataset_all->reduce(visMomCuts);
  //std::cout << "Selected " << dataset_selected->numEntries() << " TTree entries..." << std::endl;
  RooAbsData* dataset_selected = 0;

  std::map<std::string, fitManager*> fitResults;

  for ( unsigned iDecayMode = kElectron_Muon; iDecayMode <= kThreeProng1Pi0; ++iDecayMode ) {
    double iDecayMode_encoded = encodeDecayMode(iDecayMode);
    
    TString iDecayMode_string = getDecayMode_string(iDecayMode);

    TString leg1iDecayModeSign = ( iDecayMode_encoded < 0. ) ? "+" : "-";
    TString leg1iDecayModeSelection = Form("TMath::Abs(leg1DecayMode %s %2.1f) < 0.1", 
					   leg1iDecayModeSign.Data(), TMath::Abs(iDecayMode_encoded));
    if ( iDecayMode == kOneProngGt0Pi0 ) 
      leg1iDecayModeSelection = 
	Form("leg1DecayMode > %2.1f && leg1DecayMode < %2.1f", 
	     encodeDecayMode(kOneProng1Pi0) - 0.1, encodeDecayMode(kOneProng2Pi0) + 0.1);

    TString leg2iDecayModeSign = ( iDecayMode_encoded < 0. ) ? "+" : "-";
    TString leg2iDecayModeSelection = Form("TMath::Abs(leg2iDecayMode %s %2.1f) < 0.1", 
					   leg2iDecayModeSign.Data(), TMath::Abs(iDecayMode_encoded));
    if ( iDecayMode == kOneProngGt0Pi0 ) 
      leg2iDecayModeSelection = 
	Form("leg2DecayMode > %2.1f && leg2DecayMode < %2.1f", 
	     encodeDecayMode(kOneProng1Pi0) - 0.1, encodeDecayMode(kOneProng2Pi0) + 0.1);
    
    if ( iDecayMode != kElectron_Muon ) continue;

    /*****************************************************************************
     ******** Run fit for leg1 with no cuts on visible decay products applied ****
     *****************************************************************************/

    RooAbsData* dataset_leg1_all = dataset_all->reduce(leg1iDecayModeSelection.Data());

    TString fitManagerName_leg1_dR_all = Form("%s_%s_%s_%s", iDecayMode_string.Data(), "leg1", "dR", "all");
    if ( fitResults.find(fitManagerName_leg1_dR_all.Data()) == fitResults.end() ) {
      fitResults[fitManagerName_leg1_dR_all.Data()] = 
	new fitManager(&leg1Pt, &leg1VisInvisDeltaRLab, iDecayMode, 
		       Form("%s_%s_%s", "leg1", "dR", "all"), numPtBins, ptBinning);
    }
    fitManager* fitManager_leg1_dR_all = fitResults[fitManagerName_leg1_dR_all.Data()];
    TString outputFileName_leg1_dR_all = Form("plots/fitTauDecayKinePlots_%s_leg1_dR_all.eps", iDecayMode_string.Data());
    fitManager_leg1_dR_all->runPrefit(inputFileName_histograms, 
				      "VisInvisDeltaRLab", "all", outputFileName_leg1_dR_all);
    fitManager_leg1_dR_all->runFit(dataset_leg1_all, outputFileName_leg1_dR_all);
    
    TString fitManagerName_leg1_angle_all = Form("%s_%s_%s_%s", iDecayMode_string.Data(), "leg1", "angle", "all");
    if ( fitResults.find(fitManagerName_leg1_angle_all.Data()) == fitResults.end() ) {
      fitResults[fitManagerName_leg1_angle_all.Data()] = 
	new fitManager(&leg1Energy, &leg1VisInvisAngleLab, iDecayMode, 
		       Form("%s_%s_%s", "leg1", "angle", "all"), numPtBins, ptBinning);
    }
    fitManager* fitManager_leg1_angle_all = fitResults[fitManagerName_leg1_angle_all.Data()];
    TString outputFileName_leg1_angle_all = Form("plots/fitTauDecayKinePlots_%s_leg1_angle_all.eps", iDecayMode_string.Data());
    fitManager_leg1_angle_all->runPrefit(inputFileName_histograms, 
					 "VisInvisAngleLab", "all", outputFileName_leg1_angle_all);
    fitManager_leg1_angle_all->runFit(dataset_leg1_all, outputFileName_leg1_angle_all);
    
    /*****************************************************************************
     ********   Run fit for leg1 with cuts on visible decay products applied  ****
     *****************************************************************************/

    RooAbsData* dataset_leg1_selected = dataset_selected->reduce(leg1iDecayModeSelection.Data());

    TString fitManagerName_leg1_dR_selected = Form("%s_%s_%s_%s", iDecayMode_string.Data(), "leg1", "dR", "selected");
    if ( fitResults.find(fitManagerName_leg1_dR_selected.Data()) == fitResults.end() ) {
      fitResults[fitManagerName_leg1_dR_selected.Data()] = 
	new fitManager(&leg1Pt, &leg1VisInvisDeltaRLab, iDecayMode, 
		       Form("%s_%s_%s", "leg1", "dR", "selected2"), numPtBins, ptBinning);
    }
    fitManager* fitManager_leg1_dR_selected = fitResults[fitManagerName_leg1_dR_selected.Data()];
    TString outputFileName_leg1_dR_selected = Form("plots/fitTauDecayKinePlots_%s_leg1_dR_selected.eps", iDecayMode_string.Data());
    fitManager_leg1_dR_selected->runPrefit(inputFileName_histograms, 
					   "VisInvisDeltaRLab", "selected2", outputFileName_leg1_dR_selected);
    fitManager_leg1_dR_selected->runFit(dataset_leg1_selected, outputFileName_leg1_dR_selected);
    
    TString fitManagerName_leg1_angle_selected = Form("%s_%s_%s_%s", iDecayMode_string.Data(), "leg1", "angle", "selected");
    if ( fitResults.find(fitManagerName_leg1_angle_selected.Data()) == fitResults.end() ) {
      fitResults[fitManagerName_leg1_angle_selected.Data()] = 
	new fitManager(&leg1Energy, &leg1VisInvisAngleLab, iDecayMode, 
		       Form("%s_%s_%s", "leg1", "angle", "selected2"), numPtBins, ptBinning);
    }
    fitManager* fitManager_leg1_angle_selected = fitResults[fitManagerName_leg1_angle_selected.Data()];
    TString outputFileName_leg1_angle_selected = Form("plots/fitTauDecayKinePlots_%s_leg1_angle_selected.eps", iDecayMode_string.Data());
    fitManager_leg1_angle_selected->runPrefit(inputFileName_histograms, 
					      "VisInvisAngleLab", "selected2", outputFileName_leg1_angle_selected);
    fitManager_leg1_angle_selected->runFit(dataset_leg1_selected, outputFileName_leg1_angle_selected);
    
    /*****************************************************************************
     ******** Run fit for leg2 with no cuts on visible decay products applied ****
     *****************************************************************************/

    RooAbsData* dataset_leg2_all = dataset_all->reduce(leg2iDecayModeSelection.Data());

    TString fitManagerName_leg2_dR_all = Form("%s_%s_%s_%s", iDecayMode_string.Data(), "leg2", "dR", "all");
    if ( fitResults.find(fitManagerName_leg2_dR_all.Data()) == fitResults.end() ) {
      fitResults[fitManagerName_leg2_dR_all.Data()] = 
	new fitManager(&leg2Pt, &leg2VisInvisDeltaRLab, iDecayMode, 
		       Form("%s_%s_%s", "leg2", "dR", "all"), numPtBins, ptBinning);
    }
    fitManager* fitManager_leg2_dR_all = fitResults[fitManagerName_leg2_dR_all.Data()];
    TString outputFileName_leg2_dR_all = Form("plots/fitTauDecayKinePlots_%s_leg2_dR_all.eps", iDecayMode_string.Data());
    fitManager_leg2_dR_all->runPrefit(inputFileName_histograms, 
				      "VisInvisDeltaRLab", "all", outputFileName_leg2_dR_all);
    fitManager_leg2_dR_all->runFit(dataset_leg2_all, outputFileName_leg2_dR_all);
    
    TString fitManagerName_leg2_angle_all = Form("%s_%s_%s_%s", iDecayMode_string.Data(), "leg2", "angle", "all");
    if ( fitResults.find(fitManagerName_leg2_angle_all.Data()) == fitResults.end() ) {
      fitResults[fitManagerName_leg2_angle_all.Data()] = 
	new fitManager(&leg2Energy, &leg2VisInvisAngleLab, iDecayMode, 
		       Form("%s_%s_%s", "leg2", "angle", "all"), numPtBins, ptBinning);
    }
    fitManager* fitManager_leg2_angle_all = fitResults[fitManagerName_leg2_angle_all.Data()];
    TString outputFileName_leg2_angle_all = Form("plots/fitTauDecayKinePlots_%s_leg2_angle_all.eps", iDecayMode_string.Data());
    fitManager_leg2_angle_all->runPrefit(inputFileName_histograms, 
					 "VisInvisAngleLab", "all", outputFileName_leg2_angle_all);
    fitManager_leg2_angle_all->runFit(dataset_leg2_all, outputFileName_leg2_angle_all);
    
    /*****************************************************************************
     ********   Run fit for leg2 with cuts on visible decay products applied  ****
     *****************************************************************************/

    RooAbsData* dataset_leg2_selected = dataset_selected->reduce(leg2iDecayModeSelection.Data());

    TString fitManagerName_leg2_dR_selected = Form("%s_%s_%s_%s", iDecayMode_string.Data(), "leg2", "dR", "selected");
    if ( fitResults.find(fitManagerName_leg2_dR_selected.Data()) == fitResults.end() ) {
      fitResults[fitManagerName_leg2_dR_selected.Data()] = 
	new fitManager(&leg2Pt, &leg2VisInvisDeltaRLab, iDecayMode, 
		       Form("%s_%s_%s", "leg2", "dR", "selected2"), numPtBins, ptBinning);
    }
    fitManager* fitManager_leg2_dR_selected = fitResults[fitManagerName_leg2_dR_selected.Data()];
    TString outputFileName_leg2_dR_selected = Form("plots/fitTauDecayKinePlots_%s_leg2_dR_selected.eps", iDecayMode_string.Data());
    fitManager_leg2_dR_selected->runPrefit(inputFileName_histograms, 
					   "VisInvisDeltaRLab", "selected2", outputFileName_leg2_dR_selected);
    fitManager_leg2_dR_selected->runFit(dataset_leg2_selected, outputFileName_leg2_dR_selected);
    
    TString fitManagerName_leg2_angle_selected = Form("%s_%s_%s_%s", iDecayMode_string.Data(), "leg2", "angle", "selected");
    if ( fitResults.find(fitManagerName_leg2_angle_selected.Data()) == fitResults.end() ) {
      fitResults[fitManagerName_leg2_angle_selected.Data()] = 
	new fitManager(&leg2Energy, &leg2VisInvisAngleLab, iDecayMode, 
		       Form("%s_%s_%s", "leg2", "angle", "selected2"), numPtBins, ptBinning);
    }
    fitManager* fitManager_leg2_angle_selected = fitResults[fitManagerName_leg2_angle_selected.Data()];
    TString outputFileName_leg2_angle_selected = Form("plots/fitTauDecayKinePlots_%s_leg2_angle_selected.eps", iDecayMode_string.Data());
    fitManager_leg2_angle_selected->runPrefit(inputFileName_histograms, 
					      "VisInvisAngleLab", "selected2", outputFileName_leg2_angle_selected);
    fitManager_leg2_angle_selected->runFit(dataset_leg2_selected, outputFileName_leg2_angle_selected);
  }

//--- write results of prefit to ROOT file
  TFile* outputFile_prefit = new TFile("fitTauDecayKinePlots.root", "RECREATE");
  TDirectory* dirVisInvisDeltaRLab = outputFile_prefit->mkdir("VisInvisDeltaRLab");  
  TDirectory* dirVisInvisDeltaRLab_all = dirVisInvisDeltaRLab->mkdir("all");
  TDirectory* dirVisInvisDeltaRLab_selected = dirVisInvisDeltaRLab->mkdir("selected");
  TDirectory* dirVisInvisAngleLab = outputFile_prefit->mkdir("VisInvisAngleLab");
  TDirectory* dirVisInvisAngleLab_all = dirVisInvisAngleLab->mkdir("all");
  TDirectory* dirVisInvisAngleLab_selected = dirVisInvisAngleLab->mkdir("selected");
  for ( std::map<std::string, fitManager*>::iterator prefitResult = fitResults.begin();
	prefitResult != fitResults.end(); ++prefitResult ) {
    if      ( prefitResult->second->label_.Contains("_dR_")       && 
	      prefitResult->second->label_.Contains("_all_")      ) dirVisInvisDeltaRLab_all->cd();
    else if ( prefitResult->second->label_.Contains("_angle_")    && 
	      prefitResult->second->label_.Contains("_all_")      ) dirVisInvisAngleLab_all->cd();
    else if ( prefitResult->second->label_.Contains("_dR_")       && 
	      prefitResult->second->label_.Contains("_selected_") ) dirVisInvisDeltaRLab_selected->cd();
    else if ( prefitResult->second->label_.Contains("_angle_")    && 
	      prefitResult->second->label_.Contains("_selected_") ) dirVisInvisAngleLab_selected->cd();
    else assert(0);

    prefitResult->second->prefitResultSmearedLandauMP_->Write();
    prefitResult->second->prefitResultSmearedLandauWidth_->Write();
    prefitResult->second->prefitResultSmearedLandauGSigma_->Write();
    prefitResult->second->prefitResultSkewedGaussianMean_->Write();
    prefitResult->second->prefitResultSkewedGaussianSigma_->Write();
    prefitResult->second->prefitResultSkewedGaussianAlpha_->Write();
    prefitResult->second->prefitResultMix_->Write();
  }
  delete outputFile_prefit;

  ostream* outputFile_fit = new std::ofstream("fitTauDecayKinePlots.txt", std::ios::out);
  for ( std::map<std::string, fitManager*>::iterator fitResult = fitResults.begin();
	fitResult != fitResults.end(); ++fitResult ) {
    fitResult->second->writeFitResults(*outputFile_fit, fitResult->second->label_);
  }
  delete outputFile_fit;

  delete dataTree;
}
