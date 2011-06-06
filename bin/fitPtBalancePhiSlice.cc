#include <boost/program_options.hpp>

#include "TFile.h"
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
#include "TauAnalysis/CandidateTools/interface/TauDecayKinePdf.h"

using namespace RooFit;
namespace po = boost::program_options;

int main(int argc, char **argv) {
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  TCanvas canvas("blah", "blah", 600, 800);
  canvas.Divide(1, 2);

  std::string modelsFilePath;
  std::string dataFilePath;
  std::string outputFilePath;
  std::string modelToFit;
  std::string componentToPlot;
  int phiSlice;
  int legToFit;
  std::string type;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "Show help message")
    ("model-file,m",
     po::value<string>(&modelsFilePath)->default_value("models.root"),
     "model workspace file")
    ("fit,f", po::value<std::string>(&modelToFit), "Model name to fit")
    ("phiSlice,p", po::value<int>(&phiSlice)->default_value(0),
     "Index of phi slice")
    ("leg", po::value<int>(&legToFit)->default_value(1), "Leg (1 or 2) to fit")
    ("type", po::value<std::string>(&type)->default_value(""),
     "Available options: Leg1High, Leg2High, and empty")
    ("data-file,d", po::value<std::string>(&dataFilePath), "data workspace file")
    ("output-file,o", po::value<std::string>(&outputFilePath), "fit output file")
    ("component", po::value<std::string>(&componentToPlot)->default_value(""),
     "Plot sub-component of model");
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

  TFile modelsFile(modelsFilePath.c_str(), "READ");
  if (modelsFile.IsZombie()) {
    std::cerr << "Error: could not open model file " << modelsFilePath << std::endl;
    return 3;
  }
  RooWorkspace* modelWS = static_cast<RooWorkspace*>(
      modelsFile.Get("ptBalanceModeFit"));
  if (!modelWS) {
    std::cerr << "Error: could not load model workspace from "
      << modelsFilePath << std::endl;
    return 4;
  }

  TFile dataFile(dataFilePath.c_str(), "READ");
  if (dataFile.IsZombie()) {
    std::cerr << "Error: could not open data file " << dataFilePath << std::endl;
    return 3;
  }
  RooWorkspace* dataWS = static_cast<RooWorkspace*>(
      dataFile.Get("ptBalanceDataWS"));
  if (!dataWS) {
    std::cerr << "Error: could not load data workspace from "
      << dataFilePath << std::endl;
    return 4;
  }

  // Load the binned dataset to do the prefit
  std::stringstream datasetName;
  datasetName << "leg" << legToFit << type << "DataBinned";
  RooAbsData* dataset = dataWS->data(datasetName.str().c_str());
  if (!dataset) {
    std::cerr << "Couldn't find dataset with name: " << datasetName.str()
      << " in the workspace!  Dumping data workspace info:" << std::endl;
    dataWS->Print("v");
  }

  RooAbsPdf* model = modelWS->pdf(modelToFit.c_str());
  if (!model) {
    std::cerr << "Couldn't find model with name: " << modelToFit
      << " in the workspace!  Dumping model workspace info:" << std::endl;
    modelWS->Print("v");
  }

  std::stringstream observableName;
  observableName << "leg" << legToFit << "Val";

  RooRealVar* x = dataWS->var(observableName.str().c_str());
  RooRealVar* mass = dataWS->var("resonanceMass");
  RooRealVar* phi = dataWS->var("dPhi");

  assert(mass);
  assert(phi);
  assert(x);

  std::cout << "Mass variable has " << mass->getBins()
    << " bins, from " << mass->getMin() << " to " << mass->getMax()
    << std::endl;

  std::cout << "Phi variable has " << phi->getBins()
    << " bins, from " << phi->getMin() << " to " << phi->getMax()
    << std::endl;

  std::cout << " In the initial dataset there are "
    << dataset->numEntries() << " bins and "
    << dataset->sumEntries() << " entries." << std::endl;

  double phiMin = phi->getMin();
  double phiMax = phi->getMax();
  double phiBinSize = (phiMax - phiMin)/phi->getBins();
  double lowPhiEdge = phiMin + phiSlice*(phiBinSize);
  double highPhiEdge = phiMin + (phiSlice+1)*(phiBinSize);

  std::stringstream phiCut;
  phiCut << lowPhiEdge << " < dPhi && dPhi < " << highPhiEdge;

  std::cout << "Reducing to a phi slice between "
    << lowPhiEdge << " and " << highPhiEdge << std::endl;

  std::auto_ptr<RooAbsData> phiSliceData(dataset->reduce(
      RooArgSet(*x, *mass), phiCut.str().c_str()));

  std::cout << " In the phi slice dataset there are "
    << phiSliceData->numEntries() << " bins and "
    << phiSliceData->sumEntries() << " entries." << std::endl;

  std::cout << "Beginning mass loop." << std::endl;
  int nMassBins = mass->getBins();

  std::cout << "Opening summary canvas" << std::endl;
  std::stringstream summaryFileName;
  summaryFileName << "fit_"
    << modelToFit << "phi_slice_"
    << phiSlice << "_summary.ps";
  std::stringstream summaryFileNameOpen;
  summaryFileNameOpen << summaryFileName.str() << "[";
  canvas.Print(summaryFileNameOpen.str().c_str());

  std::auto_ptr<RooArgSet> components(model->getComponents());
  std::cout << "The model has " << components->getSize() << " components:"
    << std::endl;
  model->Print("v");

  TObject *comp = NULL;
  TIterator* compIter = components->createIterator();
  int iComp = 0;
  while ((comp = compIter->Next())) {
    iComp++;
    std::cout << iComp << ") " << comp->GetName() << std::endl;
    comp->Print("v");
  }

  // Check if we don't want to fit some of the parameters
  std::stringstream staticVarSetName;
  staticVarSetName << modelToFit << "_static";
  std::cout << "Checking if we want to set pars constant: "
    << staticVarSetName.str() << std::endl;
  const RooArgSet* staticVars = modelWS->set(staticVarSetName.str().c_str());
  if (staticVars) {
    std::cout << "Setting variables constant..." << std::endl;
    compIter = staticVars->createIterator();
    while ((comp = compIter->Next())) {
      RooRealVar* var = dynamic_cast<RooRealVar*>(comp);
      assert(var);
      std::cout << "Setting " << var->GetName() << " constant." << std::endl;
      var->setConstant(true);
    }
  }

  TVectorD massPoints(nMassBins);
  TVectorD massBinning(3);
  TVectorD phiEdges(2);
  phiEdges[0] = lowPhiEdge;
  phiEdges[1] = highPhiEdge;

  double massMin = mass->getMin();
  double massMax = mass->getMax();
  double massBinSize = (massMax - massMin)/mass->getBins();
  massBinning[0] = mass->getBins();
  massBinning[1] = massMin;
  massBinning[2] = massMax;

  std::cout << "Running fit!" << std::endl;
  RooFitResult* result =
    model->fitTo(*phiSliceData,
        ConditionalObservables(*mass),
        PrintEvalErrors(0), PrintLevel(-1), NumCPU(4), Save(true));

  result->Print("v");

  for (int i = 0; i < nMassBins; ++i) {
    double lowMassEdge = massMin + i*massBinSize;
    double highMassEdge = massMin + (i+1)*massBinSize;
    massPoints[i] = (highMassEdge - lowMassEdge)*0.5 + lowMassEdge;

    std::stringstream massCut;
    massCut << lowMassEdge <<
      " < resonanceMass && resonanceMass < " << highMassEdge;
    std::auto_ptr<RooAbsData> massBin(
        phiSliceData->reduce(RooArgSet(*x), massCut.str().c_str()));

    std::cout << "Plotting the " << i << "th mass slice, for mass between "
      << lowMassEdge << " and " << highMassEdge
      << ". This mass-phi region has " << massBin->sumEntries()
      << " entries." << std::endl;

    std::auto_ptr<RooArgSet> parameters(model->getParameters(*massBin));
    RooRealVar* massVar = dynamic_cast<RooRealVar*>(parameters->find("resonanceMass"));

    assert(massVar);
    massVar->setVal((highMassEdge - lowMassEdge)*0.5 + lowMassEdge);

    std::auto_ptr<RooPlot> frame(x->frame(Range(0, 2.)));
    massBin->plotOn(frame.get());
    double plotMaximum = 1.2*frame->GetMaximum();


    // Make sure we always start from the same conditions
    //modelWS->loadSnapshot("initialConditions");
    model->plotOn(frame.get(),
        //ProjWData(RooArgSet(*mass), *massBin));
        ProjWData(*massBin));

    if (componentToPlot != "") {
      model->plotOn(frame.get(), Components(componentToPlot.c_str()),
          LineColor(kRed), ProjWData(RooArgSet(*mass), *massBin));
    }
    frame->SetMaximum(plotMaximum);

    std::stringstream massRangeStr;
    massRangeStr << "(" << lowMassEdge << " - " << highMassEdge << ")";
    std::stringstream phiRangeStr;
    phiRangeStr << "(" << lowPhiEdge << " - " << highPhiEdge << ")";

    std::stringstream plotTitle;
    plotTitle << "X for mass " << massRangeStr.str()
      << " and #Delta#phi " << phiRangeStr.str();
    frame->SetTitle(plotTitle.str().c_str());

    canvas.cd(1);
    gPad->SetLogy(false);
    frame->Draw();
    canvas.cd(2);
    gPad->SetLogy(true);
    frame->Draw();
    canvas.Print(summaryFileName.str().c_str());
  }

  std::stringstream summaryFileNameClose;
  summaryFileNameClose << summaryFileName.str() << "]";
  canvas.SaveAs(summaryFileNameClose.str().c_str());
  std::cout << "Done writing summary canvas, converting to pdf..." << std::endl;
  std::stringstream converterCommand;
  converterCommand << "epstopdf " << summaryFileName.str();
  int convertResult = system(converterCommand.str().c_str());
  if (!convertResult)
    std::cout << "Converted to pdf successfully." << std::endl;
  else {
    std::cerr << "Could not convert to pdf!" << std::endl;
  }

  std::cout << "Saving output workspace" << std::endl;
  TFile outputFile(outputFilePath.c_str(), "RECREATE");
  outputFile.cd();
  RooWorkspace ws("ptBalanceModeFit");
  ws.import(*model);
  ws.import(massPoints, "massPoints");
  ws.import(massBinning, "massBinning");
  ws.import(phiEdges, "phiEdges");
  ws.Print("v");
  ws.Write();
  outputFile.Close();
  return 0;
}
