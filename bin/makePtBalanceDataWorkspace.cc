#include <boost/program_options.hpp>

#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TROOT.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"

using namespace RooFit;
namespace po = boost::program_options;

int main(int argc, char **argv) {

  // Define program options
  po::options_description desc("Allowed options");

  int nPhiSlices;
  int varPhiBins;
  int massSliceSize;
  int startMass;
  int stopMass;
  int ptBins;
  int skip;
  double maxScaledPt;
  std::string outputFile;
  std::string inputFiles;
  std::string treename;

  desc.add_options()
    ("help", "Show help message")
    ("phi-slices", po::value<int>(&nPhiSlices)->default_value(20),
     "Number of phi slices")
    ("phi-variable-bins", po::value<int>(&varPhiBins)->default_value(0),
     "If nonzero, use a variable in size for phi.  Each bin will contain "
     "an equal number of events.")
    ("mass-slice-size", po::value<int>(&massSliceSize)->default_value(10),
     "Size of mass slice in GeV")
    ("start-mass", po::value<int>(&startMass)->default_value(65),
     "Lowest mass to consider")
    ("stop-mass", po::value<int>(&stopMass)->default_value(505),
     "Highest mass to consider")
    ("max-scaledpt", po::value<double>(&maxScaledPt)->default_value(2.0),
     "Maximum range of scaled pt")
    ("pt-bins", po::value<int>(&ptBins)->default_value(200),
     "Number of bins in the observable")
    ("output-file,o", po::value<string>(&outputFile), "output file")
    ("tree-name", po::value<string>(&treename)->default_value("makePtBalanceNtuple/ptBalanceNtuple"), "TTree name")
    ("skip-selection", po::value<int>(&skip)->default_value(0), "skip kin. selection")
    ("input-files,i",
     po::value<string>(&inputFiles),
     "input file glob");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

  if (!vm.count("input-files") || !vm.count("output-file")) {
    std::cerr << "No input or output files specified! Usage:\n"
      << desc << std::endl;
    return 2;
  }

  // Load the input files.
  std::cout << "Loading data from " << inputFiles <<  std::endl;

  TChain chain(treename.c_str());
  int filesAdded = chain.Add(inputFiles.c_str());
  if (!filesAdded) {
    std::cerr << "No valid input files specified!" << std::endl;
    return 1;
  } else {
    std::cout << "Loaded " << filesAdded << " files." << std::endl;
  }

  RooRealVar leg1Pt("leg1Pt", "leg1Pt", 0, 700);
  RooRealVar leg2Pt("leg2Pt", "leg2Pt", 0, 700);
  RooRealVar leg1VisPt("leg1VisPt", "leg1VisPt", 0, 700);
  RooRealVar leg2VisPt("leg2VisPt", "leg2VisPt", 0, 700);
  RooRealVar leg1VisEta("leg1VisEta", "leg1VisEta", -5, 5);
  RooRealVar leg2VisEta("leg2VisEta", "leg2VisEta", -5, 5);
  RooRealVar leg1Leg2DPhi("leg1Leg2DPhi", "leg1Leg2DPhi", 0, 3.14159);
  RooRealVar resonanceMass("resonanceMass", "resonanceMass", 0, 1400);

  // Build the initial dataset
  RooDataSet data("data", "dataset",
      RooArgSet(
        leg1Pt, leg2Pt,
        leg1VisPt, leg2VisPt,
        leg1VisEta, leg2VisEta,
        leg1Leg2DPhi, resonanceMass),
      Import(chain));

  std::cout << "The initial dataset has " << data.numEntries() << std::endl;

  std::cout << "Applying kinematic selections..." << std::endl;

  RooDataSet* selectedData = NULL;

  if (!skip) {
    selectedData = dynamic_cast<RooDataSet*>(data.reduce(
          RooArgSet(leg1Pt, leg2Pt, leg1Leg2DPhi, resonanceMass),
          "leg1VisPt > 15 && leg2VisPt > 20 && "
          "abs(leg1VisEta) < 2.1 && abs(leg2VisEta) < 2.3"));
  } else {
    selectedData = dynamic_cast<RooDataSet*>(data.reduce(
          RooArgSet(leg1Pt, leg2Pt, leg1Leg2DPhi, resonanceMass)));
  }

  std::cout << "Kinematic selection has " << selectedData->numEntries()
    << " entries." << std::endl;

  std::cout << "Adding extra columns...";

  RooFormulaVar leg1ValFunc("leg1Val",
      "2.0*leg1Pt/resonanceMass", RooArgSet(leg1Pt, resonanceMass));

  RooFormulaVar leg2ValFunc("leg2Val",
      "2.0*leg2Pt/resonanceMass", RooArgSet(leg2Pt, resonanceMass));

  RooFormulaVar dPhiFunc("dPhi", "#pi - #Delta #phi",
      "abs(3.14159265 - leg1Leg2DPhi)", RooArgSet(leg1Leg2DPhi));

  RooRealVar* leg1Val = static_cast<RooRealVar*>(selectedData->addColumn(
        leg1ValFunc));
  RooRealVar* leg2Val = static_cast<RooRealVar*>(selectedData->addColumn(
        leg2ValFunc));
  RooRealVar* dPhi = static_cast<RooRealVar*>(selectedData->addColumn(
        dPhiFunc));

  leg1Val->setRange(0, maxScaledPt);
  leg1Val->setBins(ptBins);
  leg2Val->setRange(0, maxScaledPt);
  leg2Val->setBins(ptBins);
  dPhi->setRange(0, 2.0);
  dPhi->setBins(nPhiSlices);
  resonanceMass.setRange(startMass, stopMass);
  int nMassBins = TMath::Nint((stopMass-startMass)*1.0/massSliceSize);
  resonanceMass.setBins(nMassBins);

  std::cout << " done." << std::endl;

  std::cout << "Keeping only interesting variables..." << std::endl;
  RooAbsData* reducedSelectedData = selectedData->reduce(
      RooArgSet(*leg1Val, *leg2Val, *dPhi, resonanceMass));

  std::cout << "Creating global binned dataset" << std::endl;

  RooDataHist leg1DataBinned("leg1DataBinned", "leg1DataBinned",
      RooArgSet(*leg1Val, *dPhi, resonanceMass), *reducedSelectedData);

  RooDataHist leg2DataBinned("leg2DataBinned", "leg2DataBinned",
      RooArgSet(*leg2Val, *dPhi, resonanceMass), *reducedSelectedData);

  // Don't need it anymore.
  //delete selectedData;

  /*
  std::cout << "Splitting into leg1High and leg2High datasets" << std::endl;
  RooAbsData* leg1HighData = reducedSelectedData->reduce(
      "leg1Val > leg2Val");
  std::cout << "Leg1High has " << leg1HighData->numEntries() << std::endl;

  RooAbsData* leg2HighData = reducedSelectedData->reduce(
      "leg2Val >= leg1Val");

  std::cout << "Leg2High has " << leg2HighData->numEntries() << std::endl;

  std::cout << "Creating binned datasets for Leg 1 High" << std::endl;

  RooDataHist leg1Leg1HighDataBinned("leg1Leg1HighDataBinned",
      "leg1Leg1HighDataBinned",
      RooArgSet(*leg1Val, *dPhi, resonanceMass),
      *leg1HighData);

  RooDataHist leg2Leg1HighDataBinned("leg2Leg1HighDataBinned",
      "leg2Leg1HighDataBinned",
      RooArgSet(*leg2Val, *dPhi, resonanceMass),
      *leg1HighData);

  std::cout << "Creating binned datasets for Leg 2 High" << std::endl;

  RooDataHist leg1Leg2HighDataBinned("leg1Leg2HighDataBinned",
      "leg1Leg2HighDataBinned",
      RooArgSet(*leg1Val, *dPhi, resonanceMass),
      *leg2HighData);

  RooDataHist leg2Leg2HighDataBinned("leg2Leg2HighDataBinned",
      "leg2Leg2HighDataBinned",
      RooArgSet(*leg2Val, *dPhi, resonanceMass),
      *leg2HighData);
  */

  std::cout << "Building workspace." << std::endl;
  RooWorkspace ws("ptBalanceDataWS");

  //ws.import(*reducedSelectedData);
  ws.import(leg1DataBinned);
  ws.import(leg2DataBinned);
  //ws.import(leg1Leg1HighDataBinned);
  //ws.import(leg2Leg1HighDataBinned);
  //ws.import(leg1Leg2HighDataBinned);
  //ws.import(leg2Leg2HighDataBinned);

  ws.Print("v");

  TFile output(outputFile.c_str(), "RECREATE");
  output.cd();
  ws.Write();
  output.Close();

  return 0;
}

