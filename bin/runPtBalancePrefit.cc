#include <boost/program_options.hpp>

#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TVectorD.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"

using namespace RooFit;
namespace po = boost::program_options;

int main(int argc, char **argv) {
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  TCanvas canvas("blah", "blah", 600, 800);
  //canvas.Divide(1, 2);

  std::string prefitFilePath;
  std::string outputFilePath;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "Show help message")
    ("prefit-file",
     po::value<string>(&prefitFilePath),
     "prefit workspace file")
    ("output-file,o", po::value<std::string>(&outputFilePath), "fit output file");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

  if (!vm.count("output-file")) {
    std::cerr << "No output file specified! Usage:\n"
      << desc << std::endl;
    return 2;
  }

  TFile prefitFile(prefitFilePath.c_str(), "READ");
  if (prefitFile.IsZombie()) {
    std::cerr << "Error: could not open prefit file " << prefitFilePath << std::endl;
    return 3;
  }
  RooWorkspace* prefitWS = static_cast<RooWorkspace*>(
      prefitFile.Get("ptBalanceModeFit"));
  if (!prefitWS) {
    std::cerr << "Error: could not load model workspace from "
      << prefitFilePath << std::endl;
    return 4;
  }

  // Get the collection of fit results
  TClonesArray* fitResults = dynamic_cast<TClonesArray*>(
      prefitWS->obj("fitResults"));

  TVectorD* massPoints = dynamic_cast<TVectorD*>(
      prefitWS->obj("massPoints"));

  TVectorD* massBinning = dynamic_cast<TVectorD*>(
      prefitWS->obj("massBinning"));

  TVectorD* phiEdges = dynamic_cast<TVectorD*>(
      prefitWS->obj("phiEdges"));

  if (!fitResults) {
    std::cerr << "Error: could not load fit results from workspace."
      << std::endl;
    return 5;
  }

  std::cout << "Found " << fitResults->GetSize()
    << " mass points to fit." << std::endl;

  RooFitResult* firstFitResult = static_cast<RooFitResult*>(fitResults->First());

  assert(firstFitResult);

  //const RooArgList& pars = firstFitResult->floatParsFinal();
  RooArgList pars = firstFitResult->floatParsFinal();
  pars.add(firstFitResult->constPars());
  std::cout << "I'm going to fit the parameters: " << pars << std::endl;

  int nPars = pars.getSize();
  std::map<std::string, TH1F> parameters;

  std::stringstream summaryFileName;
  summaryFileName << "prefit_"
    << "parameter_summary.ps";
  std::stringstream summaryFileNameOpen;
  summaryFileNameOpen << summaryFileName.str() << "[";
  canvas.Print(summaryFileNameOpen.str().c_str());

  // Initialize the histograms
  for (int iPar = 0; iPar < nPars; ++iPar) {
    const RooRealVar* var = static_cast<const RooRealVar*>(pars.at(iPar));
    assert(var);
    TH1F histo(var->GetName(), var->GetName(),
        TMath::Nint((*massBinning)[0]), (*massBinning)[1], (*massBinning)[2]);

    // Now fill it with all the fit results
    int nMassPoints = fitResults->GetSize();
    for (int iMass = 0; iMass < nMassPoints; ++iMass) {
      RooFitResult* fit = static_cast<RooFitResult*>((*fitResults)[iMass]);
      assert(fit);
      RooArgList fitPars = fit->floatParsFinal();
      fitPars.add(fit->constPars());
      double mass = (*massPoints)[iMass];
      const RooRealVar* fitVar = static_cast<const RooRealVar*>(fitPars.at(iPar));
      Int_t bin = histo.FindBin(mass);
      histo.SetBinContent(bin, fitVar->getVal());
      histo.SetBinError(bin, fitVar->getError());
    }
    double minBinCont = histo.GetBinContent(histo.GetMinimumBin());
    double maxBinCont = histo.GetBinContent(histo.GetMaximumBin());
    histo.SetMaximum( (maxBinCont > 0) ? maxBinCont*1.1 : maxBinCont*0.9 );
    histo.SetMinimum( (minBinCont < 0) ? minBinCont*1.1 : minBinCont*0.9 );

    parameters[var->GetName()] = histo;
    histo.Draw();
    canvas.Print(summaryFileName.str().c_str());
  }
  std::stringstream summaryFileNameClose;
  summaryFileNameClose << summaryFileName.str() << "]";
  canvas.SaveAs(summaryFileNameClose.str().c_str());

  std::stringstream converterCommand;
  converterCommand << "epstopdf " << summaryFileName.str();
  int convertResult = system(converterCommand.str().c_str());
  if (!convertResult)
    std::cout << "Converted to pdf successfully." << std::endl;
  else {
    std::cerr << "Could not convert to pdf!" << std::endl;
  }

  std::cout << "Built histograms." << std::endl;


  return 0;
}
