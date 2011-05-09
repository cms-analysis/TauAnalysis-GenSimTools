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

#include "TauAnalysis/FittingTools/interface/RooSkewNormal.h"
#include "TauAnalysis/FittingTools/interface/TauDecayKinePdf.h"

using namespace RooFit;
namespace po = boost::program_options;

RooRealVar* makeVar(RooWorkspace& ws, const std::string& command) {
  RooRealVar* output = static_cast<RooRealVar*>(ws.factory(command.c_str()));
  return output;
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
      ws.factory("leg1Val[0, 4.0]"));
  RooRealVar* leg2Val = static_cast<RooRealVar*>(
      ws.factory("leg2Val[0, 4.0]"));
  RooRealVar* dPhi = static_cast<RooRealVar*>(
    ws.factory("dPhi[0, 2.0]"));
  RooRealVar* resonanceMass = static_cast<RooRealVar*>(
    ws.factory("resonanceMass[70, 495]"));

  // Generic location variable
  RooRealVar* location = makeVar(ws, "location[0.8, 0.5, 2]");
  RooRealVar* location2 = makeVar(ws, "location2[0.8, 0.0, 20]");
  // Generic scale variable
  RooRealVar* scale = makeVar(ws, "scale[0.5, 0, 4]");
  RooRealVar* scale2 = makeVar(ws, "scale2[0.5, 0, 4]");
  // Generic skew variable
  RooRealVar* skew = makeVar(ws, "skew[1.0, -50, 50]");
  RooRealVar* skew2 = makeVar(ws, "skew2[1.0, -50, 50]");
  skew->setError(1.0);

  // Add a Gaussian for testing
  ws.factory("RooGaussian::gaussLeg1(leg1Val, location, scale)");

  //ws.import(skewNormalLeg1);

  RooRealVar* relativeScaleFactor = makeVar(ws,
      "relativeScaleFactor[0.5, 0.2, 1.0]");

  RooRealVar* lessThanMix = makeVar(ws, "mix[0.5, 0, 0.95]");

  RooFormulaVar relativeScale("relativeScale", "relativeScale",
      "@0*@1", RooArgList(*relativeScaleFactor, *scale));

  RooSkewNormal skewNormalLeg1("skewNormalLeg1", "skewNormalLeg1",
      *leg1Val, *location, relativeScale, *skew);

  RooSkewNormal skewNormal2Leg1("skewNormal2Leg1", "skewNormal2Leg1",
      *leg1Val, *location, *scale, *skew);
  //location->setConstant(true);

  RooAddPdf doubleSkewNormalLeg1("doubleSkewNormalLeg1", "doubleSkewNormalLeg1",
      skewNormal2Leg1,
      skewNormalLeg1,
      *lessThanMix);

  RooConstVar neg2("neg2", "neg2", -3.0);
  RooExponential shrinkSkinny("shrinkSkinny", "shrinkSkinny",
      *lessThanMix, neg2);

  RooConstVar locationConstraintFactor("locConstFactor", "locConstFactor", 0.20);
  RooGaussian locationConstraint("locConstraint", "locConstraint",
      *location, *location2, locationConstraintFactor);

  RooProdPdf constrainedDoubleSkewLeg1(
      "constrainedDoubleSkewLeg1",
      "constrainedDoubleSkewLeg1",
      RooArgList(
        shrinkSkinny,
        //locationConstraint,
        doubleSkewNormalLeg1));

  RooConstVar zero("zero", "zero", 0.0);
  RooGamma gammaLeg1("gammaLeg1", "gammaLeg1",
      *leg1Val, *location2, *scale2, zero);

  RooAddPdf skewNormPlusGamma("skewNormPlusGamma", "skewNormPlusGamma",
      gammaLeg1, skewNormal2Leg1, *lessThanMix);

  RooProdPdf skewNormPlusGammaConstrained("skewNormPlusGammaConstrained",
      "skewNormPlusGammaConstrained",
      RooArgList(shrinkSkinny, skewNormPlusGamma));

  ws.import(doubleSkewNormalLeg1, RecycleConflictNodes());
  ws.import(constrainedDoubleSkewLeg1, RecycleConflictNodes());
  ws.import(skewNormPlusGammaConstrained, RecycleConflictNodes());

  // Make Leg 2
  RooSkewNormal skewNormalLeg2("skewNormalLeg2", "skewNormalLeg2",
      *leg2Val, *location, relativeScale, *skew);

  RooSkewNormal skewNormal2Leg2("skewNormal2Leg2", "skewNormal2Leg2",
      *leg2Val, *location, *scale, *skew);

  RooAddPdf doubleSkewNormalLeg2("doubleSkewNormalLeg2", "doubleSkewNormalLeg2",
      skewNormal2Leg2,
      skewNormalLeg2,
      *lessThanMix);

  RooProdPdf constrainedDoubleSkewLeg2(
      "constrainedDoubleSkewLeg2",
      "constrainedDoubleSkewLeg2",
      RooArgList(shrinkSkinny, doubleSkewNormalLeg2));

  ws.import(doubleSkewNormalLeg2, RecycleConflictNodes());
  ws.import(constrainedDoubleSkewLeg2, RecycleConflictNodes());

  // Christians PDF
  RooRealVar* gmean = makeVar(ws, "gmean[0.8, 0.1, 2]");
  RooRealVar* gsigma = makeVar(ws, "gsigma[0.1, 0.0, 2]");
  RooRealVar* alpha = makeVar(ws, "alpha[0, -20, 20]");
  alpha->setConstant(true);
  RooRealVar* offset = makeVar(ws, "offset[0, -1, 1]");
  RooRealVar* slope = makeVar(ws, "slope[2, 0, 3]");

  RooRealVar* C = makeVar(ws, "C[0.25, 0, 1]");

  RooRealVar* mp1 = makeVar(ws, "mp1[0.5, -10, 10]");
  RooRealVar* width1 = makeVar(ws, "width1[1, 0, 10]");

  RooRealVar* mp2 = makeVar(ws, "mp2[1, -10, 10]");
  mp2->setConstant(true);
  RooRealVar* width2 = makeVar(ws, "width2[1, 0, 10]");
  width2->setConstant(true);
  RooRealVar* x0 = makeVar(ws, "x0[1.1, 0.5, 10]");
  RooRealVar* dx1 = makeVar(ws, "dx1[10, 0.5, 20]");
  dx1->setConstant(true);
  x0->setConstant(true);

  TauDecayKinePdf tauDecayKinePdf("tauDecayKinePdf", "tauDecayKinePdf",
      *leg1Val, *gmean, *gsigma, *alpha, *slope, *offset, *C,
      *mp1, *width1, *mp2, *width2, *x0, *dx1);

  ws.import(tauDecayKinePdf, RecycleConflictNodes());

  // Final fits

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



  // Save the initial conditions in a snapshot.
  ws.saveSnapshot("initialConditions", ws.allVars());

  ws.Print("v");

  TFile output(outputFile.c_str(), "RECREATE");
  ws.Write();
  output.Close();
  return 0;
}
