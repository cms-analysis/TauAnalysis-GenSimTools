#include <stdio.h>
#include <iostream>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"

#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooLinkedList.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"

// PDFs
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooLandau.h"
#include "RooGExpModel.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooGamma.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"

using namespace RooFit;


class RooResonanceDecayPtPdf : public RooAbsPdf {
  public:
    RooResonanceDecayPtPdf(const char *name, const char *title,
        RooAbsReal& _scaledDecayProdPt, RooAbsReal& _maxScaledPt);
    RooResonanceDecayPtPdf(const RooResonanceDecayPtPdf& other,
        const char *name=0);
    virtual TObject* clone(const char* newname) const {
      return new RooResonanceDecayPtPdf(*this, newname);
    };
    inline virtual ~RooResonanceDecayPtPdf() {}
    Double_t evaluate() const;
    Int_t getAnalyticalIntegral(
        RooArgSet& allVars, RooArgSet& analVars, const char *rangeName) const;
    Double_t analyticalIntegral(Int_t code) const;
  private:
    // Some shortcuts
    inline Double_t square(Double_t x) const { return x*x; }
    inline Double_t sinThetaSquared() const {
      return square(decayProdSinTheta_);
    }

    inline Double_t scaledDecayProdPtSquared() const {
      return square(scaledDecayProdPt_);
    }

    inline Double_t maxScaledPtSquared() const {
      return square(maxScaledPt_);
    }

    Double_t indefiniteIntegral(Int_t code, Double_t integrandValue) const;
    RooRealProxy scaledDecayProdPt_;
    RooRealProxy maxScaledPt_;
    mutable int functionCall_;
    //ClassDef(RooResonanceDecayPtPdf, 0) // Resonance decay Pt PDF
};

RooResonanceDecayPtPdf::RooResonanceDecayPtPdf(
    const char *name, const char *title,
    RooAbsReal& _scaledDecayProdPt, RooAbsReal& _maxScaledPt):
  RooAbsPdf(name, title),
  scaledDecayProdPt_("scaledDecayProdPt",
      "Pt of decay product divided by half of the resonance mass",
      this, _scaledDecayProdPt),
  maxScaledPt_("maxScaledPt", "Scale of the decay", this,
      _maxScaledPt) {}


RooResonanceDecayPtPdf::RooResonanceDecayPtPdf(
    const RooResonanceDecayPtPdf& other, const char *name):
  RooAbsPdf(other, name),
  scaledDecayProdPt_("scaledDecayProdPt", this, other.scaledDecayProdPt_),
  maxScaledPt_("maxScaledPt", this, other.maxScaledPt_)
{}

Double_t RooResonanceDecayPtPdf::evaluateCutoffImpl(
    double pt, double maxPt) const {
  double cutoffPointVal = cutoffPoint(maxPt);
  assert(pt > cutoffPointVal);
  Double_t valueAtCutoff = evaluateNormalImpl(cutoffPointVal, maxPt);
  return valueAtCutoff*TMath::Exp(
      (cutoffPointVal - scaledDecayProdPt_)*cutOffDecay_);
}

Double_t RooResonanceDecayPtPdf::evaluateCutoffRegion() const {

}

Double_t RooResonanceDecayPtPdf::evaluateNormalImpl(
    double pt, double maxPt) const {
  Double_t xSquared = square(pt);
  Double_t sSquared = square(maxPt);
  return pt*TMath::Sqrt(sSquared - xSquared);
}

Double_t RooResonanceDecayPtPdf::evaluateNormalRegion() const {
  assert(scaledDecayProdPt_ < maxScaledPt_);
  return evaluateNormalImpl(scaledDecayProdPt_, maxScaledPt_);
}

Double_t RooResonanceDecayPtPdf::evaluate() const {
  double cutoffPointVal = cutoffPoint(maxScaledPt_);
  if (scaledDecayProdPt_ > cutoffPointVal)) {
    return evaluateCutoffRegion();
  } else
    return evaluateNormalRegion();
}

Double_t RooResonanceDecayPtPdf::indefiniteIntegralNormalRegion(Int_t code,
    Double_t integrandValue) const {
  if (code == 1) {
    assert(scaledDecayProdPt_ <= cutoffPointVal(integrandValue));
    // Integrating w.r.t. maximum scaled Pt for a fixed value of Pt
    return scaledDecayProdPt_*TMath::Log(
        2*(integrandValue + TMath::Sqrt(
            square(integrandValue) - square(scaledDecayProdPt_))));
  } else if (code == 2) {
    assert(integrandValue <= cutoffPointVal(maxScaledPt_));
    // Integrating w.r.t. scaled Pt for a fixed value of scaled Pt
    return -TMath::Sqrt(square(maxScaledPt_) - square(integrandValue));
  }
}

Double_t RooResonanceDecayPtPdf::indefiniteIntegralCutoffRegion(Int_t code,
    Double_t integrandValue) const {
  if (code == 1) {
    // Integrating w.r.t. maximum scaled Pt for a fixed value of Pt
    return scaledDecayProdPt_*TMath::Log(
        2*(integrandValue + TMath::Sqrt(
            square(integrandValue) - square(scaledDecayProdPt_))));
  } else if (code == 2) {
    assert(integrandValue <= cutoffPointVal(maxScaledPt_));
    // Integrating w.r.t. scaled Pt for a fixed value of scaled Pt
    return -TMath::Sqrt(square(maxScaledPt_) - square(integrandValue));
  }
}

Int_t RooResonanceDecayPtPdf::getAnalyticalIntegral(
    RooArgSet& allVars, RooArgSet& analVars, const char *rangeName) const {
  //if (matchArgs(allVars, analVars, scaledDecayProdPt_, maxScaledPt_)) {
    //return 3;
  //}
  if (matchArgs(allVars, analVars, scaledDecayProdPt_)) {
    return 2;
  }
  if (matchArgs(allVars, analVars, maxScaledPt_)) {
    return 1;
  }
  return 0;
}


Double_t RooResonanceDecayPtPdf::analyticalIntegral(Int_t code) const {
  if (code == 1) {
    //Integrate w.r.t scaled decay Pt
    Double_t max = scaledDecayProdPt_.max();
    Double_t maxPossible = TMath::Sqrt(sinThetaSquared() + betaGammaSquared());
    max = (max > maxPossible) ? maxPossible : max;
    Double_t min = scaledDecayProdPt_.min();
    min = min >= 0 ? min : 0;
    min = min > maxPossible ? maxPossible : min;
    return indefiniteIntegral(1, max) - indefiniteIntegral(1, min);
  }
  if (code == 2) {
    Double_t min = decayProdSinTheta_.min();
    Double_t max = decayProdSinTheta_.max();
    assert(min >= -1);
    assert(max <= 1);
    return indefiniteIntegral(2, max) - indefiniteIntegral(2, min);
  }
  if (code == 3) {
    Double_t min = betaGamma_.min();
    Double_t max = betaGamma_.max();
    assert(min >= 0);
    return indefiniteIntegral(3, max) - indefiniteIntegral(3, min);
  }
  assert(false);
  return 0;
}

int main(int argc, const char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: ./run_fit file_glob" << std::endl;
    return 1;
  }
  gROOT->SetBatch(true);
  // Load the data
  std::cout << "Loading data from " << argv[1] <<  std::endl;
  TChain chain("makePtBalanceNtuple/ptBalanceNtuple");
  chain.Add(argv[1]);

  TTree* dataTree = static_cast<TTree*>(&chain);

  /************************************************************************
   ********       Define variabes and data                             ****
   ************************************************************************/

  RooRealVar leg1Pt("leg1Pt", "leg1Pt", 0, 250);
  RooRealVar leg2Pt("leg2Pt", "leg2Pt", 0, 250);
  RooRealVar leg1VisPt("leg1VisPt", "leg1VisPt", 0, 200);
  RooRealVar leg2VisPt("leg2VisPt", "leg2VisPt", 0, 200);
  RooRealVar leg1VisEta("leg1VisEta", "leg1VisEta", -5, 5);
  RooRealVar leg2VisEta("leg2VisEta", "leg2VisEta", -5, 5);
  RooRealVar leg1Leg2DPhi("leg1Leg2DPhi", "leg1Leg2DPhi", 0, 3.14159);
  RooRealVar resonanceMass("resonanceMass", "resonanceMass", 0, 500);
  RooRealVar resonancePt("resonancePt", "resonancePt", 0, 500);

  RooFormulaVar betaTFormula("betaT",
      "sqrt(resonancePt*resonancePt/(resonanceMass*resonanceMass + resonancePt*resonancePt))",
      RooArgSet(resonanceMass, resonancePt));

  RooFormulaVar gammaTFormula("gammaT",
      "1.0/(sqrt(1 - betaT*betaT))-1.0", betaTFormula);

  RooFormulaVar betaGammaTFormula("betaGammaT",
      "betaT*gammaT", RooArgSet(betaTFormula, gammaTFormula));

  RooFormulaVar leg1VisTheta("leg1VisTheta",
      "2*atan(exp(-leg1VisEta))", leg1VisEta);

  RooFormulaVar leg1VisSinThetaFunc("leg1VisSinTheta",
      "sin(leg1VisTheta)", leg1VisTheta);


  RooFormulaVar scaledLeg1PtFunc("scaledLeg1Pt",
      "2.0*leg1Pt/resonanceMass", RooArgSet(leg1Pt, resonanceMass));

  TCanvas canvas("blah", "blah", 800, 600);

  RooDataSet data("data", "dataset",
      RooArgSet(
        leg1Pt, leg2Pt,
        leg1VisPt, leg2VisPt,
        leg1VisEta, leg2VisEta,
        leg1Leg2DPhi, resonanceMass, resonancePt),
      Import(*dataTree));

  RooRealVar* betaGammaT = static_cast<RooRealVar*>(data.addColumn(betaGammaTFormula));
  betaGammaT->setRange(0, 100);
  RooRealVar* scaledLeg1Pt = static_cast<RooRealVar*>(
      data.addColumn(scaledLeg1PtFunc));
  scaledLeg1Pt->setRange(0.005, 20);
  RooRealVar* leg1VisSinTheta = static_cast<RooRealVar*>(
      data.addColumn(leg1VisSinThetaFunc));
  leg1VisSinTheta->setRange(-1, 1);

  RooFormulaVar leg1SinThetaSqPlusBetaGSqFunc("leg1SinThetaSqPlusBetaGSqFunc",
      "sqrt(leg1VisSinTheta*leg1VisSinTheta + betaGammaT*betaGammaT)",
      RooArgSet(*leg1VisSinTheta,*betaGammaT));

  RooRealVar* leg1SinThetaSqPlusBetaGSq = static_cast<RooRealVar*>(
    data.addColumn(leg1SinThetaSqPlusBetaGSqFunc));

  RooRealVar* gammaT = static_cast<RooRealVar*>(
      data.addColumn(gammaTFormula));
  gammaT->setRange(0, 20);

  std::cout << "Initial data set has " << data.numEntries() << std::endl;

  RooAbsData* dataAnaSelection = data.reduce(
      "leg1VisPt > 15 && leg2VisPt > 20"
      "&& abs(leg1VisEta) < 2.1 && abs(leg2VisEta) < 2.3");

  std::cout << "After analysis selection data has "
    << dataAnaSelection->numEntries() << "entries" <<  std::endl;

  RooAbsData* dataBackToBack = dataAnaSelection->reduce(
      "leg1Leg2DPhi > 2.8");

  RooAbsData* dataNotBackToBack = dataAnaSelection->reduce(
      "leg1Leg2DPhi < 2.2 && leg1Leg2DPhi > 2.0");

  std::cout << "There are " << dataBackToBack->numEntries()
    << " back-to-back entries" << std::endl;

  RooRealVar betaGammaFit("betaGammaFit", "testBetaGammaFit", 0.5, 0, 10);
  RooRealVar sinThetaFit("testST", "testST", 0.98, 0, 1);
  //sinThetaFit.setConstant(true);

  RooRealVar betaGammaMean("betaGammaMean", "betaGammaMean", 0.5, 0, 30);
  RooRealVar betaGammaSigma("betaGammaSigma", "betaGammaSigma", 0.5, 0, 30);

  RooLandau betaGammaPdf("betaGammaPdf", "betaGammaPdf",
      betaGammaFit, betaGammaMean, betaGammaSigma);

  RooResonanceDecayPtPdf mypdf("mypdf", "mypdf",
      *scaledLeg1Pt, *leg1VisSinTheta, *betaGammaT);

  RooRealVar testThing("blahblah", "blahblah", 0.05, 0, 5);

  RooGaussian testicles("blah", "blah",
      *scaledLeg1Pt, *leg1SinThetaSqPlusBetaGSq, testThing);


  RooRealVar gammaDecayBtB("gammaFitBtB", "gammaFitBtB", -1, -100, 100);
  RooRealVar gammaDecayNotBtB("gammaFitNotBtB", "gammaFitNotBtB", 5, -1000, 1000);
  RooRealVar gammaDecayExpTerm("gammaFitExpTerm", "gammaFitExpTerm", 1, -1000, 1000);
  RooRealVar gammaDecayExpTerm2("gammaFitExpTerm2", "gammaFitExpTerm2", 1, 0.5, 4);

  RooGenericPdf myTestFitasdf("test", "test", "pow(@0, @2)*exp(@1*pow(@0, @3))",
      RooArgList(*gammaT, gammaDecayNotBtB, gammaDecayExpTerm, gammaDecayExpTerm2));

  RooLandau myTestFit("test2", "test2",
      *gammaT, gammaDecayNotBtB, gammaDecayExpTerm);

  RooGExpModel myTestFitFuck("test2", "test2",
      *gammaT, gammaDecayNotBtB, gammaDecayExpTerm);

  myTestFit.fitTo(*dataNotBackToBack);
  //myTestFit.fitTo(*dataBackToBack);
  //gammaFitBtB.fitTo(*dataBackToBack);

  //RooFFTConvPdf model("model", "model", betaGammaFit,
      //mypdf, betaGammaPdf);

  //model.fitTo(*testSlice,
      //ConditionalObservables(*leg1VisSinTheta)
      //);

  std::cout << "Painting graphs" << std::endl;

  RooPlot* gammaTFrame = gammaT->frame(Range(0, 2.));
  dataNotBackToBack->plotOn(gammaTFrame);
  //dataBackToBack->plotOn(gammaTFrame);
  //gammaFitBtB.plotOn(gammaTFrame);
  myTestFit.plotOn(gammaTFrame);
  gammaTFrame->Draw();
  canvas.SaveAs("gammaT.png");

  RooPlot* sinThetaFrame = leg1VisSinTheta->frame(Range(-2, 2.));
  dataNotBackToBack->plotOn(sinThetaFrame);
  sinThetaFrame->Draw();
  canvas.SaveAs("sinTheta.png");

  RooPlot* betaGammaFrame = betaGammaT->frame(Range(-0, 5.));
  dataNotBackToBack->plotOn(betaGammaFrame);
  betaGammaFrame->Draw();
  canvas.SaveAs("betaGamma.png");

  RooPlot* scaledLeg1Frame = scaledLeg1Pt->frame(Range(0, 3.0));
  //testicles.plotOn(scaledLeg1Frame,
      //ProjWData(RooArgSet(*leg1SinThetaSqPlusBetaGSq), *dataNotBackToBack, false));
  *leg1VisSinTheta = 1.0;
  *betaGammaT = 1.0;
  //mypdf.plotOn(scaledLeg1Frame,
      //ProjWData(RooArgSet(*leg1VisSinTheta, *betaGammaT), *dataNotBackToBack, true));
  dataNotBackToBack->plotOn(scaledLeg1Frame);
  mypdf.plotOn(scaledLeg1Frame, ProjWData(RooArgSet(*leg1VisSinTheta, *betaGammaT), *dataNotBackToBack), Normalization(dataNotBackToBack->numEntries()));
  scaledLeg1Frame->Draw();
  scaledLeg1Frame->Print("v");
  canvas.SaveAs("scaledLeg1Pt.png");
  return 0;
}
