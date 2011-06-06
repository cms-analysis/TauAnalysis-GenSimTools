
#include "RooHist.h"
#include "RooCurve.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TList.h>
#include <TH1.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TPaveText.h>

#include <vector>

void duplicate(TObjArray& labels, Int_t idx)
{
  assert(idx < labels.GetEntries());
  labels.Add(labels.At(idx));
}

double getIntegral(TGraph* graph)
{
  double retVal = 0.;

  Int_t numPoints = graph->GetN();
  for ( Int_t iPoint = 0; iPoint < numPoints; ++iPoint ) {
    double x, y;
    graph->GetPoint(iPoint, x, y);
    retVal += y;
  }
   
  return retVal;
}

double getMaximum(TGraph* graph)
{
  double retVal = -1.;
  
  Int_t numPoints = graph->GetN();
  for ( Int_t iPoint = 0; iPoint < numPoints; ++iPoint ) {
    double x, y;
    graph->GetPoint(iPoint, x, y);
    if ( y > retVal ) retVal = y;
  }
   
  return retVal;
}

void scale(TGraph* graph, double scaleFactor)
{
  Int_t numPoints = graph->GetN();
  for ( Int_t iPoint = 0; iPoint < numPoints; ++iPoint ) {
    double x, y;
    graph->GetPoint(iPoint, x, y);
    graph->SetPoint(iPoint, x, y*scaleFactor);

    if ( dynamic_cast<TGraphAsymmErrors*>(graph) ) {
      TGraphAsymmErrors* graphAsymmErrors = dynamic_cast<TGraphAsymmErrors*>(graph);
      double yErrUp   = graphAsymmErrors->GetErrorYhigh(iPoint);
      double yErrDown = graphAsymmErrors->GetErrorYlow(iPoint);
      graphAsymmErrors->SetPointEYhigh(iPoint, yErrUp*scaleFactor);
      graphAsymmErrors->SetPointEYlow(iPoint, yErrDown*scaleFactor);
    } else if ( dynamic_cast<TGraphErrors*>(graph) ) {
      TGraphErrors* graphErrors = dynamic_cast<TGraphErrors*>(graph);
      double xErr = graphErrors->GetErrorX(iPoint);
      double yErr = graphErrors->GetErrorY(iPoint);
      graphErrors->SetPointError(iPoint, xErr, yErr*scaleFactor);
    } 
  }
}

void niceTauDecayKinePlots()
{
  TString inputFilePath = "/afs/cern.ch/user/v/veelken/scratch0/CMSSW_4_1_3/src/TauAnalysis/GenSimTools/bin/plots/";

  TObjArray inputFileNames;
  inputFileNames.Add(new TObjString("fitTauDecayKinePlots_Electron_Muon_leg1_angle_all_Energy45to46_log.root"));
  inputFileNames.Add(new TObjString("fitTauDecayKinePlots_OneProng0Pi0_leg2_angle_all_Energy46to47_log.root"));
  inputFileNames.Add(new TObjString("fitTauDecayKinePlots_OneProngGt0Pi0_leg2_angle_all_Energy45to46_log.root"));
  inputFileNames.Add(new TObjString("fitTauDecayKinePlots_ThreeProng0Pi0_leg2_angle_all_Energy45to46_log.root"));
  inputFileNames.Add(new TObjString("fitTauDecayKinePlots_Electron_Muon_leg1_angle_selected_Energy45to50_log.root"));
  inputFileNames.Add(new TObjString("fitTauDecayKinePlots_OneProng0Pi0_leg2_angle_selected_Energy45to50_log.root"));
  inputFileNames.Add(new TObjString("fitTauDecayKinePlots_OneProngGt0Pi0_leg2_angle_selected_Energy45to50_log.root"));
  inputFileNames.Add(new TObjString("fitTauDecayKinePlots_ThreeProng0Pi0_leg2_angle_selected_Energy45to50_log.root"));
  
  TObjArray labelsDecayMode;
  labelsDecayMode.Add(new TObjString("#tau^{-} #rightarrow #mu^{-} #bar{#nu}_{#mu} #nu_{#tau}, #tau^{-} #rightarrow e^{-} #bar{#nu}_{e} #nu_{#tau}"));
  labelsDecayMode.Add(new TObjString("#tau^{-} #rightarrow #pi^{-} #nu_{#tau}"));
  labelsDecayMode.Add(new TObjString("#tau^{-} #rightarrow #pi^{-} #pi^{0} #nu_{#tau}, #tau^{-} #rightarrow #pi^{-} #pi^{0} #pi^{0} #nu_{#tau}"));
  labelsDecayMode.Add(new TObjString("#tau^{-} #rightarrow #pi^{-} #pi^{+} #pi^{-} #nu_{#tau}"));
  duplicate(labelsDecayMode, 0);
  duplicate(labelsDecayMode, 1);
  duplicate(labelsDecayMode, 2);
  duplicate(labelsDecayMode, 3);
    
  TObjArray labelsSelection;
  labelsSelection.Add(new TObjString("before p_{T}^{vis} > 15 GeV cut"));
  labelsSelection.Add(new TObjString("before p_{T}^{vis} > 20 GeV cut"));
  duplicate(labelsSelection, 1);
  duplicate(labelsSelection, 1);
  labelsSelection.Add(new TObjString("after p_{T}^{vis} > 15 GeV cut"));
  labelsSelection.Add(new TObjString("after p_{T}^{vis} > 20 GeV cut"));
  duplicate(labelsSelection, 5);
  duplicate(labelsSelection, 5);
  
  std::vector<double> xOffsetsLabel;
  xOffsetsLabel.push_back(0.45);
  xOffsetsLabel.push_back(0.45);
  xOffsetsLabel.push_back(0.19);
  xOffsetsLabel.push_back(0.19);
  xOffsetsLabel.push_back(0.45);
  xOffsetsLabel.push_back(0.39);
  xOffsetsLabel.push_back(0.19);
  xOffsetsLabel.push_back(0.19);

  std::vector<double> yOffsetsLabel;
  yOffsetsLabel.push_back(0.51);
  yOffsetsLabel.push_back(0.51);
  yOffsetsLabel.push_back(0.22);
  yOffsetsLabel.push_back(0.22);
  yOffsetsLabel.push_back(0.51);
  yOffsetsLabel.push_back(0.22);
  yOffsetsLabel.push_back(0.22);
  yOffsetsLabel.push_back(0.22);

  std::cout << "#inputFileNames  = " << inputFileNames.GetEntries()  << std::endl;
  std::cout << "#labelsDecayMode = " << labelsDecayMode.GetEntries() << std::endl;
  std::cout << "#labelsSelection = " << labelsSelection.GetEntries() << std::endl;
  std::cout << "#xOffsetsLabel   = " << xOffsetsLabel.size()        << std::endl;
  std::cout << "#yOffsetsLabel   = " << yOffsetsLabel.size()        << std::endl;

  assert(inputFileNames.GetEntries() == labelsDecayMode.GetEntries());
  assert(inputFileNames.GetEntries() == labelsSelection.GetEntries());
  assert(inputFileNames.GetEntries() == (int)xOffsetsLabel.size());
  assert(inputFileNames.GetEntries() == (int)yOffsetsLabel.size());

  unsigned numInputFiles = inputFileNames.GetEntries();
  for ( unsigned iInputFile = 0; iInputFile < numInputFiles; ++iInputFile ) {
    TObjString* inputFileName = (TObjString*)inputFileNames.At(iInputFile);

    TFile* inputFile = new TFile(TString(inputFilePath).Append(inputFileName->GetString()).Data());

    TCanvas* canvas = (TCanvas*)inputFile->Get("canvas");

    TList* objects = canvas->GetListOfPrimitives();
    
    TGraph*  histogram   = 0;
    TGraph*  fitFunction = 0;
    TLegend* legend      = 0;
    
    TIter next(objects);
    TObject* object = 0;
    unsigned iObject = 0;
    while ( (object = next()) ) {
      std::cout << "object #" << iObject << ": ClassName = " << object->ClassName() << std::endl;
      if ( std::string("RooHist")  == object->ClassName() ) histogram   = (TGraph*)object;
      if ( std::string("RooCurve") == object->ClassName() ) fitFunction = (TGraph*)object;
      if ( std::string("TLegend")  == object->ClassName() ) legend      = (TLegend*)object;
      ++iObject;
    }
    
    std::cout << "histogram   = " << histogram   << std::endl;
    std::cout << "fitFunction = " << fitFunction << std::endl;
    
    double scaleFactor = 1./getIntegral(histogram);
    scale(histogram,   scaleFactor);
    scale(fitFunction, scaleFactor);

    TCanvas* new_canvas = new TCanvas("canvas", "canvas", 1, 1, 800, 600);
    new_canvas->SetFillColor(10);
    new_canvas->SetBorderSize(2);

    new_canvas->SetLeftMargin(0.14);
    new_canvas->SetBottomMargin(0.14);

    new_canvas->SetLogy(true);

    TH1* dummyHistogram = new TH1F("dummyHistogram", "dummyHistogram", 50, 0., 25.);
    dummyHistogram->SetTitle("");
    dummyHistogram->SetStats(false);

    dummyHistogram->SetMinimum(0.5*scaleFactor);
    dummyHistogram->SetMaximum(2.*TMath::Max(getMaximum(histogram), getMaximum(fitFunction)));    

    TAxis* xAxis = dummyHistogram->GetXaxis();
    xAxis->SetTitle("#alpha(p_{vis}, p_{mis}) * E_{#tau} [GeV * rad]");
    xAxis->SetTitleOffset(1.3);
    
    TAxis* yAxis = dummyHistogram->GetYaxis();
    yAxis->SetTitle("a.u.");
    yAxis->SetTitleOffset(1.4);

    dummyHistogram->Draw();

    histogram->Draw("psame");
    fitFunction->Draw("lsame");

    legend->Draw();

    double xOffsetLabel = xOffsetsLabel[iInputFile];
    double yOffsetLabel = yOffsetsLabel[iInputFile];

    TObjString* labelDecayMode = (TObjString*)labelsDecayMode.At(iInputFile);
    TPaveText* labelDecayMode_pave = 
      new TPaveText(xOffsetLabel, yOffsetLabel + 0.07, xOffsetLabel + 0.44, yOffsetLabel + 0.13, "brNDC");
    labelDecayMode_pave->AddText(labelDecayMode->GetString().Data());
    labelDecayMode_pave->SetBorderSize(0);
    labelDecayMode_pave->SetFillColor(10);
    labelDecayMode_pave->Draw();
  
    TObjString* labelSelection = (TObjString*)labelsSelection.At(iInputFile);
    TPaveText* labelSelection_pave = 
      new TPaveText(xOffsetLabel, yOffsetLabel + 0.00, xOffsetLabel + 0.44, yOffsetLabel + 0.06, "brNDC");
    labelSelection_pave->AddText(labelSelection->GetString().Data());
    labelSelection_pave->SetBorderSize(0);
    labelSelection_pave->SetFillColor(10);
    labelSelection_pave->Draw();

    TString outputFileName = inputFileName->GetString();
    outputFileName.ReplaceAll(".root", "_nice.pdf");

    new_canvas->SaveAs(outputFileName.Data());

    delete new_canvas;

    delete inputFile;
  }
}
