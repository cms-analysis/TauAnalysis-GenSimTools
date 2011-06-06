#include <boost/program_options.hpp>

#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TVectorD.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "TauAnalysis/CandidateTools/interface/RooAbsEstimatablePdf.h"
#include "TauAnalysis/CandidateTools/interface/RooCruijff.h"

using namespace RooFit;
namespace po = boost::program_options;

namespace {
  TGraphErrors makeXYGraph(const RooDataSet& data,
      RooAbsArg& x, RooAbsArg& y) {
    std::cout << "<makeXYGraph>" << std::endl;
    double sumYError = 0;
    double xError = -1;
    double ySum = 0;
    size_t numEntries = data.numEntries();
    TGraphErrors output(numEntries);
    for(size_t i = 0; i < numEntries; ++i) {
      const RooArgSet* row = data.get(i);
      RooRealVar* xVal = dynamic_cast<RooRealVar*>(row->find(x.GetName()));
      RooRealVar* yVal = dynamic_cast<RooRealVar*>(row->find(y.GetName()));
      assert(xVal);
      assert(yVal);

      output.SetPoint(i, xVal->getVal(), yVal->getVal());
      std::cout << "Adding x: " << xVal->getVal() << " y: " << yVal->getVal()
        << std::endl;
      sumYError += yVal->getError();
      xError = xVal->getError(); // always the same
      ySum += yVal->getVal();
    }
    // Make sure we don't have bins with tiny errors driving the fit
    for(size_t i = 0; i < numEntries; ++i) {
      output.SetPointError(i, xError,
          TMath::Max(sumYError/numEntries, TMath::Abs(ySum/numEntries)));
    }
    return output;
  }
}

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
    ("models-file,m",
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
      modelsFile.Get("ptBalanceModelWS"));
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
  x->setRange(0, 2.0);
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
  summaryFileName << "prefit_"
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
  }

  // Define a dataset
  RooArgSet modelObservables(*x);

  std::auto_ptr<RooArgSet> modelParameters(model->getParameters(
      &modelObservables));

  // Add the resonance mass variable to our data
  std::auto_ptr<RooArgSet> fitResultVariables(new RooArgSet(*modelParameters));

  fitResultVariables->add(*mass);

  std::cout << "Found the following dependent variables:" << std::endl;
  fitResultVariables->Print("v");

  std::cout << "Building fit results data store..." << std::endl;

  RooDataSet fitResults("fitResult", "fitResults",
      *fitResultVariables, StoreError(*fitResultVariables));

  std::cout << " Data store has the structure: " << std::endl;
  fitResults.Print("v");

  TClonesArray results("RooFitResult", nMassBins);
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

  for (int i = 0; i < nMassBins; ++i) {
    double lowMassEdge = massMin + i*massBinSize;
    double highMassEdge = massMin + (i+1)*massBinSize;
    massPoints[i] = (highMassEdge - lowMassEdge)*0.5 + lowMassEdge;

    std::stringstream massCut;
    massCut << lowMassEdge <<
      " < resonanceMass && resonanceMass < " << highMassEdge;
    std::auto_ptr<RooAbsData> massBin(
        phiSliceData->reduce(RooArgSet(*x), massCut.str().c_str()));

    std::cout << "Fitting the " << i << "th mass slice, for mass between "
      << lowMassEdge << " and " << highMassEdge
      << ". This mass-phi region has " << massBin->sumEntries()
      << " entries." << std::endl;

    std::auto_ptr<RooPlot> frame(x->frame(Range(0, 2.)));
    massBin->plotOn(frame.get());
    double plotMaximum = 1.2*frame->GetMaximum();

    // Make sure we always start from the same conditions
    modelWS->loadSnapshot("initialConditions");
    // Try to estimate the paramters
    std::cout << "Trying to esimate parameters..." << std::endl;
    RooAbsEstimatablePdf* estModel = dynamic_cast<RooAbsEstimatablePdf*>(model);
    if (estModel) {
      std::cout << "Estimating..." << std::endl;
      RooArgSet pars = estModel->estimateParameters(*massBin);
      std::cout << "Estimated values for: " << pars << std::endl;
      std::cout << "Plotting pre-estimation" << std::endl;
      // Plot pre-estimation
      model->plotOn(frame.get(), LineStyle(kDashed));
    } else {
      std::cout << "PDF does not support pre-estimation" << std::endl;
    }


    std::auto_ptr<RooArgSet> observables(model->getObservables(*massBin));
    RooRealVar* modelObservable = dynamic_cast<RooRealVar*>(
        observables->first());
    modelObservable->setMax(2.0);

    RooFitResult* result  = model->fitTo(
        *massBin, Save(true), PrintLevel(-1), PrintEvalErrors(0));

    result->Print("v");
    std::cout << "Covariance quality: " << result->covQual() << std::endl;
    results[i] = static_cast<TObject*>(result);
    model->plotOn(frame.get());
    if (componentToPlot != "") {
      model->plotOn(frame.get(), Components(componentToPlot.c_str()),
          LineColor(kRed));
    }
    frame->SetMaximum(plotMaximum);

    if (result->covQual() == 3) {
      RooArgSet dataThisFit;
      *mass = (highMassEdge - lowMassEdge)*0.5 + lowMassEdge;
      mass->setError(massBinSize*0.5);
      dataThisFit.add(*mass);

      const RooArgList& floatParsFinal = result->floatParsFinal();
      const RooArgList& constPars = result->constPars();
      for (int iPar = 0; iPar < floatParsFinal.getSize(); ++iPar) {
        RooRealVar* par = dynamic_cast<RooRealVar*>(floatParsFinal.at(iPar));
        assert(par);
        //par->setError(TMath::Abs(par->getVal()));
        dataThisFit.add(*par);
      }
      for (int iPar = 0; iPar < constPars.getSize(); ++iPar) {
        RooRealVar* par = dynamic_cast<RooRealVar*>(constPars.at(iPar));
        assert(par);
        //par->setError(0.1*TMath::Abs(par->getVal()));
        dataThisFit.add(*par);
      }
      dataThisFit.Print("v");
      fitResults.add(dataThisFit);
    }

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

  fitResults.Print("v");

  //canvas.cd(1);
  //gPad->SetLogy(false);
  //canvas.cd(2);
  //gPad->SetLogy(false);
  //canvas.cd(0);

  std::stringstream extendModelName;
  extendModelName << modelToFit << "_fun";
  std::cout << "Finding extended function: "
    << extendModelName.str() << std::endl;

  RooAbsPdf* extendedModel = modelWS->pdf(extendModelName.str().c_str());
  if (extendedModel) {
    std::cout << "Dumping model!" << std::endl;
    extendedModel->printTree(std::cout);
  } else {
    std::cout << "Couldn't find the model!" << std::endl;
  }


  std::cout << "Printing fit trends..." << std::endl;
  TObject *modelPar = NULL;
  TIterator* modelParIter = modelParameters->createIterator();
  int iPar = 0;
  std::vector<RooPlot*> garbage;
  while ((modelPar = modelParIter->Next())) {
    std::auto_ptr<RooPlot> massFrame(mass->frame());
    RooRealVar* yvar = static_cast<RooRealVar*>(modelPar);

    std::cout << "Making trends for: " <<  yvar->GetName() << std::endl;

    fitResults.plotOnXY(massFrame.get(), YVar(*yvar));
    massFrame->SetTitle(yvar->GetName());
    //garbage.push_back(massFrame);

    std::cout << "Checking for paramterization function" << std::endl;
    std::stringstream funName;
    funName << yvar->GetName() << "_fun";
    RooAbsReal* function = modelWS->function(funName.str().c_str());
    if (function) {
      std::cout << "Found function!" << std::endl;
      fitResults.Print("v");
      function->Print("v");
      std::cout << "observables: " <<
        *(function->getObservables(fitResults)) << std::endl;
      std::cout << "parameters: " <<
        *(function->getParameters(fitResults)) << std::endl;

      std::auto_ptr<RooArgSet> func_observables(
          function->getObservables(fitResults));
      std::auto_ptr<RooArgSet> func_parameters(
          function->getParameters(fitResults));

      RooArgList func_observables_List(*func_observables);
      RooArgList func_parameters_List(*func_parameters);

      std::cout << "Building TF1 " << std::endl;
      std::auto_ptr<TF1> functionToFit(function->asTF(
          func_observables_List,
          func_parameters_List));

      std::cout << "Making TGraph" << std::endl;
      TGraphErrors graph = makeXYGraph(fitResults, *mass, *yvar);

      function->plotOn(massFrame.get(), LineColor(kRed));

      std::cout << "Fitting TGraph" << std::endl;
      graph.Fit(functionToFit.get());

      std::cout << "Copying parameter errors" << std::endl;
      for (int iFittedPar = 0; iFittedPar < functionToFit->GetNpar(); ++iFittedPar) {
        RooRealVar* var = dynamic_cast<RooRealVar*>(func_parameters_List.at(iFittedPar));
        assert(var);
        std::cout << "Par " << iFittedPar << ") "
          << functionToFit->GetParName(iFittedPar)
          // sanity check
          << ":" << var->GetName()
          << " = " << functionToFit->GetParameter(iFittedPar)
          << " +- " << functionToFit->GetParError(iFittedPar) << std::endl;

        var->setVal(functionToFit->GetParameter(iFittedPar));
        var->setError(functionToFit->GetParError(iFittedPar));
        var->Print("v");
      }

      function->plotOn(massFrame.get());
    }

    canvas.cd(1);
    massFrame->Draw();

    double minimum = massFrame->GetMinimum();
    std::cout << "Minimum is: " <<  minimum << std::endl;

    std::auto_ptr<TH1> histo(fitResults.createHistogram(yvar->GetName()));
    int firstNonZeroBin = -1;
    int lastNonZeroBin = -1;
    for (int iBin = 1; iBin < histo->GetNbinsX(); ++iBin) {
      if (histo->GetBinContent(iBin) > 0 && firstNonZeroBin != -1) {
        firstNonZeroBin = iBin;
      }
      if (histo->GetBinContent(iBin) > 0) {
        lastNonZeroBin = iBin;
      }
    }

    double firstNonZeroX = histo->GetBinCenter(firstNonZeroBin);
    double lastNonZeroX = histo->GetBinCenter(lastNonZeroBin);
    massFrame->SetMinimum(firstNonZeroX - TMath::Abs(0.1*firstNonZeroX));
    massFrame->SetMaximum(lastNonZeroX + TMath::Abs(0.1*lastNonZeroX));

    if (massFrame->GetMinimum() > 0.0) {
      canvas.cd(2);
      massFrame->Draw();
      canvas.Print(summaryFileName.str().c_str());
    } else {
      // Else we need to add an offset to make it positive everywhere
      /*
      std::stringstream formula;
      formula << "@0 - " << massFrame->GetMinimum();
      RooFormulaVar addOffset("addOffset", "addOffset",
          formula.str().c_str(), RooArgList(*yvar));
      RooRealVar* offsetVar = static_cast<RooRealVar*>(
          fitResults.addColumn(addOffset));
      std::auto_ptr<RooPlot> massFrameForOffset(mass->frame());
      fitResults.plotOn(massFrameForOffset.get(), YVar(*offsetVar));
      canvas.cd(2);
      massFrameForOffset->Draw();
      */
      canvas.Print(summaryFileName.str().c_str());
    }

    iPar++;
  }

  if (extendedModel) {
    std::cout << "Dumping model again!" << std::endl;
    extendedModel->printTree(std::cout);
  } else {
    std::cout << "Couldn't find the model!" << std::endl;
  }

  std::stringstream summaryFileNameClose;
  summaryFileNameClose << summaryFileName.str() << "]";
  canvas.SaveAs(summaryFileNameClose.str().c_str());
  std::cout << "Done writing summary canvas, converting to pdf..." << std::endl;
  std::stringstream converterCommand;
  converterCommand << "epstopdf " << summaryFileName.str();
  std::stringstream removeCommand;
  removeCommand << "rm " << summaryFileName.str();
  int convertResult = system(converterCommand.str().c_str());
  if (!convertResult) {
    std::cout << "Converted to pdf successfully, deleting .ps file." << std::endl;
    system(removeCommand.str().c_str());
  } else {
    std::cerr << "Could not convert to pdf!" << std::endl;
  }

  std::cout << "Saving output workspace" << std::endl;

  TFile outputFile(outputFilePath.c_str(), "RECREATE");
  outputFile.cd();
  RooWorkspace ws("ptBalanceModeFit");
  ws.import(*model);
  ws.import(results, "fitResults");
  ws.import(massPoints, "massPoints");
  ws.import(massBinning, "massBinning");
  ws.import(phiEdges, "phiEdges");
  if (extendedModel)
    ws.import(*extendedModel);

  ws.Print("v");
  ws.Write();
  outputFile.Close();
  return 0;
}
