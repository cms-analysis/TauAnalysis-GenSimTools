#include <boost/program_options.hpp>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

/* Command line example:
 *
 * Filter all root files in current directory, adding suffix _selected.  The
 * grep [0-5] is to reduce the total size.
 *
 * ls *root | grep -v -e "selected.root" | sed "s|\(.*\).root|-i \1.root -o \1_selected.root|" | grep -e "_[0-5]_" | xargs -P 4 -n 4 filterPtBalanceNtuples
 */

namespace po = boost::program_options;

int main(int argc, char **argv) {

  // Define program options
  po::options_description desc("Allowed options");

  std::string outputFilePath;
  std::string inputFiles;
  std::string outputTreeName;
  std::string treeName;
  std::string cut;

  desc.add_options()
    ("help", "Show help message")
    ("output-file,o", po::value<std::string>(&outputFilePath), "output file")
    ("input-files,i", po::value<std::string>(&inputFiles), "input file glob")
    ("tree,t",
     po::value<std::string>(&treeName)->default_value(
       "makePtBalanceNtuple/ptBalanceNtuple"), "input tree")
    ("outputTree",
     po::value<std::string>(&outputTreeName)->default_value("ptBalanceNtuple"),
     "output tree name")
    ("cut",
     po::value<std::string>(&cut)->default_value(
      "leg1VisPt > 15 && leg2VisPt > 20 && "
      "abs(leg1VisEta) < 2.1 && abs(leg2VisEta) < 2.3"), "selection cut");

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

  TChain inputChain(treeName.c_str());

  int filesAdded = inputChain.Add(inputFiles.c_str());
  if (!filesAdded) {
    std::cerr << "No valid input files specified!" << std::endl;
    return 1;
  } else {
    std::cout << "Loaded " << filesAdded << " files." << std::endl;
  }

  std::cout << "Input data has " << inputChain.GetEntries() << " entries."
    << std::endl;

  TFile outputFile(outputFilePath.c_str(), "RECREATE");
  outputFile.cd();

  std::cout << "Beginning selection..." << std::endl;
  TTree* copy = inputChain.CopyTree(cut.c_str());

  assert(copy);

  std::cout << "Selected data has " << copy->GetEntries() << " entries."
    << std::endl;

  copy->SetName(outputTreeName.c_str());

  std::cout << "Writing selected tree to " << outputFilePath << std::endl;
  copy->Write();

  outputFile.Close();

  return 0;
}
