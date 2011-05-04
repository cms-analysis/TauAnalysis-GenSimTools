
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

#include <iostream>
#include <iomanip>

int main(int argc, const char* argv[])
{
  if ( argc < 6 ) {
    std::cerr << "Usage: ./testRooWorkSpace_loading inputFileName wsName pdfName momName sepName" << std::endl;
    return 1;
  }

  std::cout << "<testRooWorkSpace_loading>:" << std::endl;

  TString inputFileName = argv[1];
  std::cout << " inputFileName = " << inputFileName.Data() << std::endl;
  TString wsName = argv[2];
  std::cout << " wsName = " << wsName.Data() << std::endl;
  TString pdfName = argv[3];
  std::cout << " pdfName = " << pdfName.Data() << std::endl;
  TString momName = argv[4];
  std::cout << " momName = " << momName.Data() << std::endl;
  TString sepName = argv[5];
  std::cout << " sepName = " << sepName.Data() << std::endl;
  TString outputFileName = Form("testRooWorkSpace_loading_%s", inputFileName.Data());
  outputFileName = outputFileName.ReplaceAll(".root", ".png");
  std::cout << " outputFileName = " << outputFileName.Data() << std::endl;
  
  TFile* inputFile = new TFile(inputFileName.Data());
  std::cout << "--> inputFile = " << inputFile << std::endl;

  RooWorkspace* ws = (RooWorkspace*)inputFile->Get(wsName.Data());
  std::cout << "--> ws = " << ws << std::endl;
  //ws->Print();

  RooAbsPdf* model = ws->pdf(pdfName.Data());
  std::cout << "--> model = " << model << std::endl;

  RooRealVar* mom = ws->var(momName.Data());
  std::cout << "--> mom = " << mom << std::endl;

  RooRealVar* sep = ws->var(sepName.Data());
  std::cout << "--> sep = " << sep << std::endl;

  std::vector<double> momTestValues;
  momTestValues.push_back(15.5);
  momTestValues.push_back(25.5);
  momTestValues.push_back(45.5);
  momTestValues.push_back(75.5);
  momTestValues.push_back(125.5);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 1, 1, 600, 900);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  for ( std::vector<double>::const_iterator momTestValue = momTestValues.begin();
	momTestValue != momTestValues.end(); ++momTestValue ) {
    mom->setVal(*momTestValue);

    RooDataSet* toyData = model->generate(*sep, 10000); 

    TString frameTitle = Form("P_{T} = %2.1f GeV", *momTestValue);
    RooPlot* frame = sep->frame(RooFit::Title(frameTitle.Data()), RooFit::Bins(120));

    toyData->plotOn(frame);
      
    model->plotOn(frame, RooFit::ProjWData(RooArgSet(*sep), *toyData));
    
    frame->Draw();

    canvas->SetLogy(false);
    canvas->Update();
    TString outputFileName_linear = outputFileName;
    outputFileName_linear = outputFileName_linear.ReplaceAll(".png", Form("Pt%2.1fGeV_linear.png", *momTestValue));
    outputFileName_linear = outputFileName_linear.ReplaceAll(".5GeV", "_5GeV");
    canvas->SaveAs(outputFileName_linear.Data());

    canvas->SetLogy(true);
    canvas->Update();
    TString outputFileName_log = outputFileName;
    outputFileName_log = outputFileName_linear.ReplaceAll(".png", Form("Pt%2.1fGeV_log.png", *momTestValue));
    outputFileName_log = outputFileName_log.ReplaceAll(".5GeV", "_5GeV");
    canvas->SaveAs(outputFileName_log.Data());
  }
      
  delete inputFile;

  delete canvas;
}
