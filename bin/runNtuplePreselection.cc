
#include "TauAnalysis/GenSimTools/bin/tauDecayKineAuxFunctions.h"

#include <TString.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TROOT.h>

#include <iostream>
#include <iomanip>

int main(int argc, const char* argv[])
{
  if ( argc < 5 ) {
    std::cerr << "Usage: ./runNtuplePreselection inputFileNames treeName selection outputFileName" << std::endl;
    return 1;
  }

  printTimeStamp("<runNtuplePreselection::main (begin)>");

  TString inputFileNames = argv[1];
  TString treeName = argv[2];
  TString selection = argv[3];
  TString outputFileName = argv[4];
  
  gROOT->SetBatch(true);

  // Load the data
  std::cout << "Loading data from " << inputFileNames <<  std::endl;
  TChain* inputTree = new TChain(treeName.Data());
  inputTree->Add(inputFileNames);
  std::cout << " entries = " << inputTree->GetEntries() << std::endl;

  // Copy entries passing (pre)selection
  TTree* outputTree = inputTree->CopyTree(selection.Data());
  std::cout << " selection = " << selection.Data() << std::endl;
  std::cout << "--> entries passing selection = " << outputTree->GetEntries() << std::endl;

  // Save selected data 
  TFile* outputFile = new TFile(outputFileName.Data(), "RECREATE");
  outputTree->Write();
  delete outputFile;

  delete inputTree;

  printTimeStamp("<runNtuplePreselection::main (end)>");
}
