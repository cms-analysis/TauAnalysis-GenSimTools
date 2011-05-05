
#include "TauAnalysis/GenSimTools/bin/tauDecayKineAuxFunctions.h"
#include "TauAnalysis/FittingTools/interface/TauDecayKinePdf.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFit.h"
#include "RooWorkspace.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"

#include <TString.h>
#include <TFile.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLegend.h>

#include <iostream>
#include <iomanip>

int main(int argc, const char* argv[])
{
  if ( argc < 6 ) {
    std::cerr << "Usage: ./testRooWorkSpace_loading inputFileName wsName pdfName momName sepName" << std::endl;
    return 1;
  }

  std::cout << "<testRooWorkSpace_loading>:" << std::endl;

  TString inputFileName_ws = argv[1];
  std::cout << " inputFileName_ws = " << inputFileName_ws.Data() << std::endl;
  TString wsName = argv[2];
  std::cout << " wsName = " << wsName.Data() << std::endl;
  TString pdfName = argv[3];
  std::cout << " pdfName = " << pdfName.Data() << std::endl;
  TString momName = argv[4];
  std::cout << " momName = " << momName.Data() << std::endl;
  TString sepName = argv[5];
  std::cout << " sepName = " << sepName.Data() << std::endl;
  TString inputFileName_ws_base = inputFileName_ws;
  if ( inputFileName_ws_base.Last('/') != kNPOS ) inputFileName_ws_base.Replace(0, inputFileName_ws_base.Last('/') + 1, "", 0);
  TString outputFileName = Form("testRooWorkSpace_loading_%s", inputFileName_ws_base.Data());
  outputFileName = outputFileName.ReplaceAll(".root", ".png");
  std::cout << " outputFileName = " << outputFileName.Data() << std::endl;
  
  TFile* inputFile_ws = new TFile(inputFileName_ws.Data());
  std::cout << "--> inputFile_ws = " << inputFile_ws << std::endl;

  RooWorkspace* ws = (RooWorkspace*)inputFile_ws->Get(wsName.Data());
  std::cout << "--> ws = " << ws << std::endl;
  //ws->Print();

  RooAbsPdf* model = ws->pdf(pdfName.Data());
  std::cout << "--> model = " << model << std::endl;

  RooRealVar* mom = ws->var(momName.Data());
  std::cout << "--> mom = " << mom << std::endl;

  RooRealVar* sep = ws->var(sepName.Data());
  std::cout << "--> sep = " << sep << std::endl;

  TString inputFileName_histograms = "makeTauDecayKinePlots.root";
  
  TString inputDirName = "";
  if      ( inputFileName_ws.Contains("_dR_")    ) inputDirName = "VisInvisDeltaRLabTimesPt";
  else if ( inputFileName_ws.Contains("_angle_") ) inputDirName = "VisInvisAngleLabTimesEnergy";
  assert(!inputDirName.IsNull());

  TString decayMode = "";
  for ( int iDecayMode = kElectron_Muon; iDecayMode <= kThreeProng1Pi0; ++iDecayMode ) {
    TString iDecayMode_string = getDecayMode_string(iDecayMode);
    if ( inputFileName_ws.Contains(iDecayMode_string.Data()) ) decayMode = iDecayMode_string;
  }
  assert(!decayMode.IsNull());

  TString label = "";
  if      ( inputFileName_ws.Contains("_all_")      ) label = "all";
  else if ( inputFileName_ws.Contains("_selected_") ) label = "selected";
  assert(!label.IsNull());

  TString inputDirName_full = Form("%s/%s/%s", inputDirName.Data(), label.Data(), decayMode.Data());
  std::cout << " inputDirName_full = " << inputDirName_full.Data() << std::endl;

  TFile* inputFile_histograms = new TFile(inputFileName_histograms.Data());

  std::vector<double> momTestValues;
  momTestValues.push_back(15.5);
  momTestValues.push_back(25.5);
  momTestValues.push_back(45.5);
  momTestValues.push_back(75.5);
  momTestValues.push_back(125.5);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 1, 1, 800, 600);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  canvas->SetLeftMargin(0.14);

  unsigned numToys = 10000;

  for ( std::vector<double>::const_iterator momTestValue = momTestValues.begin();
	momTestValue != momTestValues.end(); ++momTestValue ) {
    mom->setVal(*momTestValue);

    RooDataSet* toyData = model->generate(*sep, numToys); 
    toyData->SetName("toyData");

    TString frameTitle = Form("P_{T} = %2.1f GeV", *momTestValue);
    RooPlot* frame = sep->frame(RooFit::Title(frameTitle.Data()), RooFit::Bins(120));

    TString xAxisTitle = "";
    if      ( inputDirName.Contains("DeltaR") ) xAxisTitle = "P_{T}^{#tau} #cdot #DeltaR [GeV]";
    else if ( inputDirName.Contains("Angle")  ) xAxisTitle = "E_{#tau} #cdot Angle [GeV #cdot radians]";
    assert(!xAxisTitle.IsNull());

    TAxis* xAxis = frame->GetXaxis();
    xAxis->SetTitle(xAxisTitle.Data());
    xAxis->SetTitleOffset(1.2);

    TAxis* yAxis = frame->GetYaxis();
    yAxis->SetTitle("Entries");
    yAxis->SetTitleOffset(1.6);

    toyData->plotOn(frame, RooFit::MarkerStyle(21), RooFit::MarkerColor(2), RooFit::LineColor(2));

    model->plotOn(frame, RooFit::ProjWData(RooArgSet(*sep), *toyData));
    
    frame->Draw();

    TString momName = "";
    if      ( inputDirName.Contains("DeltaR") ) momName = "Pt";
    else if ( inputDirName.Contains("Angle")  ) momName = "Energy";
    assert(!momName.IsNull());

    double momMin = TMath::Floor(*momTestValue);
    double momMax = TMath::Ceil(*momTestValue);
    TString histogramName = 
      Form("%s_%s%s%2.0fto%2.0f", decayMode.Data(), label.Data(), momName.Data(), momMin, momMax);
    TH1* histogram = (TH1*)inputFile_histograms->Get(TString(inputDirName_full.Data()).Append("/").Append(histogramName.Data()));      
    std::cout << " histogramName = " << histogramName.Data() << ": histogram = " << histogram << std::endl;

    histogram->Scale(numToys/histogram->Integral());
    histogram->SetMarkerStyle(20);
    histogram->SetMarkerColor(1);
    histogram->SetLineColor(1);
    histogram->Draw("epsame");

    TLegend legend(0.68, 0.69, 0.88, 0.88, "", "brNDC");
    legend.SetBorderSize(0);
    legend.SetFillColor(0);
    TObject* graphModel   = 0;
    TObject* graphToyData = 0;
    for ( int iItem = 0; iItem < frame->numItems(); ++iItem ) {
      TString itemName = frame->nameOf(iItem);
      TObject* item = frame->findObject(itemName.Data());
      if      ( itemName.Contains("pdf")     ) graphModel   = item;
      else if ( itemName.Contains("toyData") ) graphToyData = item; 
    }
    if ( graphToyData != 0 ) legend.AddEntry(graphToyData, "Toy Data", "p");
    if ( graphModel   != 0 ) legend.AddEntry(graphModel,   "Fit",      "l");
    legend.AddEntry(histogram, "TAUOLA", "p");
    legend.Draw();

    canvas->SetLogy(false);
    canvas->Update();
    TString outputFileName_linear = outputFileName;
    outputFileName_linear = outputFileName_linear.ReplaceAll(".png", Form("Pt%2.1fGeV_linear.png", *momTestValue));
    outputFileName_linear = outputFileName_linear.ReplaceAll(".5GeV", "_5GeV");
    canvas->SaveAs(outputFileName_linear.Data());

    canvas->SetLogy(true);
    canvas->Update();
    TString outputFileName_log = outputFileName;
    outputFileName_log = outputFileName_log.ReplaceAll(".png", Form("Pt%2.1fGeV_log.png", *momTestValue));
    outputFileName_log = outputFileName_log.ReplaceAll(".5GeV", "_5GeV");
    canvas->SaveAs(outputFileName_log.Data());
  }
  
  delete inputFile_ws;
  delete inputFile_histograms;

  delete canvas;
}
