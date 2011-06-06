#include <boost/program_options.hpp>

#include "TFile.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooGamma.h"
#include "RooConstVar.h"

#include "TauAnalysis/CandidateTools/interface/RooSkewNormal.h"
#include "TauAnalysis/CandidateTools/interface/TauDecayKinePdf.h"
#include "TauAnalysis/CandidateTools/interface/RooCruijff.h"

using namespace RooFit;
namespace po = boost::program_options;

RooRealVar* makeVar(RooWorkspace& ws, const std::string& command) {
  RooRealVar* output = static_cast<RooRealVar*>(ws.factory(command.c_str()));
  return output;
}

std::string makeErf(const std::string& name, const std::string& depvar,
    double xoffset, double valAtOffset, double xPlateau, double valAtPlateu) {
  double scaleFactor = (valAtPlateu - valAtOffset);
  double turnon = 1.0/(xPlateau - xoffset);
  std::stringstream output;
  output << "expr::" << name << "("
    << "'@1*TMath::Erf( (@0-@2 -" << xoffset << ")*@3 ) + @4',"
    << "{" << depvar << ","
    << name << "_yscale[" << scaleFactor << ", -20, 20],"
    << name << "_xshift[0, -40, 40],"
    << name << "_turnon[" << turnon << ", 0, 1],"
    << name << "_yoffset[" << valAtOffset << ", -10, 10]})";
  std::cout << "makeErf: " << output.str() << std::endl;
  return output.str();
}


std::string makePolynomial(const std::string& name, const std::string& depvar,
    double x0, double y0, double x1, double y1, double min=-10, double max=10) {
  double slope_term = (y1 - y0)/(x1 - x0);
  double offset_term = y1 - slope_term*x1;

  std::stringstream output;
  output << "RooPolyVar::" << name << "(" << depvar << ", "
    << "{"
    << name << "_a[" << offset_term << ", " << min << ", " << max << "], "
    << name << "_b[" << slope_term << ", " << min << ", " << max << "]})";
  std::cout << "PolyVar: " << output.str() << std::endl;
  return output.str();
}

int main(int argc, char **argv) {
  std::string outputFile;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "Show help message")
    ("output-file,o", po::value<string>(&outputFile), "output file");
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

  RooWorkspace ws("ptBalanceModelWS");

  // Build the dependent variables
  RooRealVar* leg1Val = static_cast<RooRealVar*>(
      ws.factory("leg1Val[0, 2.0]"));
  RooRealVar* leg2Val = static_cast<RooRealVar*>(
      ws.factory("leg2Val[0, 2.0]"));
  RooRealVar* dPhi = static_cast<RooRealVar*>(
    ws.factory("dPhi[0, 2.0]"));
  RooRealVar* resonanceMass = static_cast<RooRealVar*>(
    ws.factory("resonanceMass[70, 495]"));

  ws.factory("RooSkewNormal::sn_leg1("
    "leg1Val, sn_location[0.8, 0.5, 2], sn_scale[0.5, 0, 4],"
    "sn_skew[1.0, -50, 50])");

  ws.factory(makeErf("sn_location_fun", "resonanceMass",
        60, 0.78, 250, 0.64).c_str());

  ws.factory(makeErf("sn_scale_fun", "resonanceMass",
        60, 0.3, 120, 0.4).c_str());

  ws.factory(makeErf("sn_skew_fun", "resonanceMass",
        60, 2.5, 450, 1.7).c_str());

  ws.factory("RooSkewNormal::sn_leg1_fun("
    "leg1Val, sn_location_fun, sn_scale_fun,"
    "sn_skew_fun)");

  ws.factory("TauDecayKinePdf::tdk_leg1("
    "leg1Val,"
    "tdk_gmean[1, 0.1, 3],"
    "tdk_gsigma[0.1, 0.0, 2],"
    //"tdk_alpha[0, -20, 20],"
    "tdk_alpha[0, -20, 0],"
    "tdk_slope[0],"
    "tdk_offset[0],"
    "tdk_C[1.0],"
    "tdk_mp1[0.5, 0, 10],"
    "tdk_width1[1, 1e-3, 10],"
    "tdk_mp2[0.5, 0, 10],"
    "tdk_width2[1, 1e-3, 10],"
    "tdk_x0[1.1, 0.5, 10],"
    "tdk_dx1[0.5, 0.05, 1])"
    );


  std::cout << "Building gmean function" << std::endl;

  //ws.factory(makePolynomial("tdk_gmean_fun", "resonanceMass",
        //100, 1, 450, 1.8).c_str());
  ws.factory(makeErf("tdk_gmean_fun", "resonanceMass",
        60, 1.1, 500, 1.58).c_str());

  std::cout << "Building gsigma function" << std::endl;
  //ws.factory(makePolynomial("tdk_gsigma_fun", "resonanceMass",
  //      100, 0.2, 450, 0.4).c_str());
  ws.factory(makeErf("tdk_gsigma_fun", "resonanceMass",
        60, 0.19, 500, 0.38).c_str());

  std::cout << "Building alpha function" << std::endl;
  //ws.factory(makePolynomial("tdk_alpha_fun", "resonanceMass",
  //      100, -0.65, 450, -0.35).c_str());
  ws.factory(makeErf("tdk_alpha_fun", "resonanceMass",
        60, -0.65, 225, -0.4).c_str());

  /*
  std::cout << "Building slope function" << std::endl;
  ws.factory(makePolynomial("tdk_slope_fun", "resonanceMass",
        100, 0, 450, 0).c_str());

  std::cout << "Building offset function" << std::endl;
  ws.factory(makePolynomial("tdk_offset_fun", "resonanceMass",
        100, 0, 450, 0).c_str());

  std::cout << "Building C function" << std::endl;
  ws.factory(makePolynomial("tdk_C_fun", "resonanceMass",
        100, 1, 450, 1).c_str());
        */

  std::cout << "Building mp1 function" << std::endl;
  ws.factory(makePolynomial("tdk_mp1_fun", "resonanceMass",
        100, 0.99, 450, 1, -1, 2).c_str());

  ws.factory(makePolynomial("tdk_width1_fun", "resonanceMass",
        100, 1e-4, 450, 1e-4).c_str());

  ws.factory(makePolynomial("tdk_mp2_fun", "resonanceMass",
        100, 1.04, 450, 1.02).c_str());

  ws.factory(makePolynomial("tdk_width2_fun", "resonanceMass",
        100, 0.03, 450, 0.01).c_str());

  ws.factory(makePolynomial("tdk_x0_fun", "resonanceMass",
        100, 1.02, 450, 1.0).c_str());

  ws.factory(makePolynomial("tdk_dx1_fun", "resonanceMass",
        100, 0.16, 450, 0.16, 0, 0.5).c_str());

  ws.factory("TauDecayKinePdf::tdk_leg1_fun("
    "leg1Val,"
    "tdk_gmean_fun,"
    "tdk_gsigma_fun,"
    //"tdk_alpha[0, -20, 20],"
    "tdk_alpha_fun,"
    "tdk_slope[0],"
    "tdk_offset[0],"
    "tdk_C[1.0],"
    "tdk_mp1_fun,"
    "tdk_width1_fun,"
    "tdk_mp2_fun,"
    "tdk_width2_fun,"
    "tdk_x0_fun,"
    "tdk_dx1_fun)"
    );

  //ws.factory("RooCruijff::cru_leg1(leg1Val,"
      //"cru_mean[1, 0, 2], cru_sigmaL[0.1, 0, 2], cru_sigmaR[0.1, 0, 2],"
      //"cru_alphaL[0.1, 0, 2], cru_alphaR[0.1, 0, 2])");

  // Define parameters that shouldn't be included in the final fit
  std::cout << "Defining constant parameters" << std::endl;
  ws.defineSet("tdk_leg1_fun_static", "tdk_x0_fun_a");
  ws.extendSet("tdk_leg1_fun_static", "tdk_x0_fun_b");
  ws.extendSet("tdk_leg1_fun_static", "tdk_dx1_fun_a");
  ws.extendSet("tdk_leg1_fun_static", "tdk_dx1_fun_b");



  /*
  // Location of skew normal
  std::cout << "Building Skew Normal location" << std::endl;
  ws.factory(
      "RooPolyVar::sn_loc_fun(resonanceMass,"
      "{loc_a[1, 0.5, 2], loc_b[0, -0.01, 0.01]})");

  // Scale of skew normal
  std::cout << "Building SN scale function" << std::endl;
  ws.factory(
      "RooPolyVar::sn_scale_fun(resonanceMass,"
      "{scale_a[0.15, 0.0, 2], scale_b[0, -0.01, 0.01]})");

  // Skew of skew normal (minimum @ -40)
  std::cout << "Building SN skew function" << std::endl;
  ws.factory(
      "expr::sn_skew_fun("
      "'0.5*(@1+40.0)*(1 - tanh(@2*(@0 - @3))) - 40.0',"
      "resonanceMass, skew_start[1, -5, 5], skew_drop[0.25, 0, 1], skew_trans[450, 100, 600])");

  std::cout << "Building gamma gamma paramter function" << std::endl;
  ws.factory(
      "RooPolyVar::gamma_gamma_fun(resonanceMass,"
      "{gamma_a[7, 1, 40], gamma_b[0.0125, -0.05, 0.05]})");

  std::cout << "Building gamma beta function" << std::endl;
  ws.factory(
      "RooPolyVar::gamma_beta_fun(resonanceMass,"
      "{beta_a[0.1, 0.01, 2], beta_b[0, -0.01, 0.01]})");

  std::cout << "Building gamma mix function" << std::endl;
  ws.factory(
      "expr::sn_gamma_mix_fun("
      "'0.5*(@1-@4)*(1 - tanh(@2*(@0 - @3))) + @4',"
      "resonanceMass, mix_start[0.6, 0.00, 1], mix_turnoff[0.01, 0, 1.0],"
      "mix_turnloc[325, 100, 500], mix_end[0.1, 0, 1])");

  std::cout << "Building full fit" << std::endl;
  ws.factory(
      "SUM::full_fit("
      "sn_gamma_mix_fun*"
      "RooGamma::gamma(leg1Val, gamma_gamma_fun, gamma_beta_fun, 0.0),"
      "RooSkewNormal::sn(leg1Val, sn_loc_fun, sn_scale_fun, sn_skew_fun))");
      */

  // Save the initial conditions in a snapshot.
  ws.saveSnapshot("initialConditions", ws.allVars());

  ws.Print("v");

  TFile output(outputFile.c_str(), "RECREATE");
  ws.Write();
  output.Close();
  return 0;
}
