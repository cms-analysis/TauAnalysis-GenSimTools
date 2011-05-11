
#include "TauAnalysis/GenSimTools/bin/tauDecayKineAuxFunctions.h"

#include <TString.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TROOT.h>
#include <TMath.h>
#include <TH1.h>

#include <iostream>
#include <iomanip>

TH1* bookHistogram(const TString& decayMode, const TString& variable, const TString& label, double momMin, double momMax)
{
  TString momName_string  = "Mom";
  TString momTitle_string = "Momentum";
  if        ( variable.Contains("DeltaR") ) {
    momName_string  = "Pt";
    momTitle_string = "P_{T}";
  } else if ( variable.Contains("Angle")  ) {
    momName_string  = "Energy";
    momTitle_string = "Energy";
  }

  TString histogramName = 
    Form("%s_%s%s%2.0fto%2.0f", decayMode.Data(), label.Data(), momName_string.Data(), momMin, momMax);
  TString histogramTitle = 
    Form("%s %s %s: %s%2.0f..%2.0f", variable.Data(), decayMode.Data(), label.Data(), momTitle_string.Data(), momMin, momMax);
  bool isTimesMom = (variable.Contains("TimesPt") || variable.Contains("TimesEnergy"));
  double xMax = ( isTimesMom ) ? 12.0 : 1.0;
  TH1* histogram = new TH1F(histogramName.Data(), histogramTitle.Data(), 120, 0., xMax);
  return histogram;
}

struct plotEntryOneMomBin
{
  plotEntryOneMomBin(const TString& variable, const TString& label, double momMin, double momMax)
    : momMin_(momMin),
      momMax_(momMax)
  {
    sepElectron_Muon_  = bookHistogram(getDecayMode_string(kElectron_Muon),  variable, label, momMin, momMax);
    sepOneProng0Pi0_   = bookHistogram(getDecayMode_string(kOneProng0Pi0),   variable, label, momMin, momMax);
    sepOneProng1Pi0_   = bookHistogram(getDecayMode_string(kOneProng1Pi0),   variable, label, momMin, momMax);
    sepOneProng2Pi0_   = bookHistogram(getDecayMode_string(kOneProng2Pi0),   variable, label, momMin, momMax);
    sepOneProngGt0Pi0_ = bookHistogram(getDecayMode_string(kOneProngGt0Pi0), variable, label, momMin, momMax);
    sepThreeProng0Pi0_ = bookHistogram(getDecayMode_string(kThreeProng0Pi0), variable, label, momMin, momMax);
    sepThreeProng1Pi0_ = bookHistogram(getDecayMode_string(kThreeProng1Pi0), variable, label, momMin, momMax);
  }

  ~plotEntryOneMomBin()
  {
    delete sepElectron_Muon_;
    delete sepOneProng0Pi0_;
    delete sepOneProng1Pi0_;
    delete sepOneProng2Pi0_;
    delete sepOneProngGt0Pi0_;
    delete sepThreeProng0Pi0_;
    delete sepThreeProng1Pi0_;
  }

  void fill(double mom, double sep, int decayMode)
  {
    if ( mom > momMin_ && mom < momMax_ ) {
      TH1* histogram = getHistogram(decayMode);
      if ( histogram ) histogram->Fill(sep);
      if ( decayMode == kOneProng1Pi0 || decayMode == kOneProng2Pi0 ) sepOneProngGt0Pi0_->Fill(sep);
    }
  }

  void Draw(int decayMode, const char* option, int color = 1)
  {
    TH1* histogram = getHistogram(decayMode);
    if ( histogram ) {
      histogram->SetStats(false);      
      if ( !histogram->GetSumw2N() ) histogram->Sumw2();
      histogram->Scale(1./histogram->Integral());
      histogram->SetMaximum(1.);
      histogram->SetMarkerColor(color);
      //histogram->SetMarkerStyle(20);
      //histogram->SetMarkerSize(1.0);
      histogram->SetLineColor(color);
      //histogram->SetLineWidth(2);
      histogram->Draw(option);
    }
  }

  TH1* getHistogram(int decayMode)
  {
    if      ( decayMode == kElectron_Muon  ) return sepElectron_Muon_;
    else if ( decayMode == kOneProng0Pi0   ) return sepOneProng0Pi0_;
    else if ( decayMode == kOneProng1Pi0   ) return sepOneProng1Pi0_;
    else if ( decayMode == kOneProng2Pi0   ) return sepOneProng2Pi0_;
    else if ( decayMode == kOneProngGt0Pi0 ) return sepOneProngGt0Pi0_;
    else if ( decayMode == kThreeProng0Pi0 ) return sepThreeProng0Pi0_;
    else if ( decayMode == kThreeProng1Pi0 ) return sepThreeProng1Pi0_;
    else return 0;
  }

  double momMin_;
  double momMax_;

  TH1* sepElectron_Muon_;
  TH1* sepOneProng0Pi0_;
  TH1* sepOneProng1Pi0_;
  TH1* sepOneProng2Pi0_;
  TH1* sepOneProngGt0Pi0_;
  TH1* sepThreeProng0Pi0_;
  TH1* sepThreeProng1Pi0_;
};

struct plotEntryAllMomBins
{
  plotEntryAllMomBins(const TString& variable, const TString& label, const TArrayD& momBinning, const TArrayD& sepBinning)
  {
    int numMomBins = momBinning.GetSize() - 1;
    for ( int iMomBin = 0; iMomBin < numMomBins; ++iMomBin ) {
      double momMin = momBinning[iMomBin];
      double momMax = momBinning[iMomBin + 1];
      plotEntryOneMomBin* newPlotEntry = new plotEntryOneMomBin(variable, label, momMin, momMax);
      plotEntries_.push_back(newPlotEntry);
    }
  }
  
  ~plotEntryAllMomBins()
  {
    for ( std::vector<plotEntryOneMomBin*>::iterator it = plotEntries_.begin();
	  it != plotEntries_.end(); ++it ) {
      delete (*it);
    }
  }

  void fill(double mom, double sep, int genDecayMode)
  {
    for ( std::vector<plotEntryOneMomBin*>::iterator plotEntry = plotEntries_.begin();
	  plotEntry != plotEntries_.end(); ++plotEntry ) {
      (*plotEntry)->fill(mom, sep, genDecayMode);
    }
  }

  plotEntryOneMomBin* findPlotEntry(double mom)
  {
    plotEntryOneMomBin* retVal = 0;

    for ( std::vector<plotEntryOneMomBin*>::iterator plotEntry = plotEntries_.begin();
	  plotEntry != plotEntries_.end(); ++plotEntry ) {
      if ( (*plotEntry)->momMin_ <= mom && (*plotEntry)->momMax_ >= mom ) {
	retVal = (*plotEntry);
	break;
      }
    }

    return retVal;
  }

  void Draw(TCanvas* canvas, const TString& outputFileName, const char* option)
  {
    for ( unsigned iDecayMode = kElectron_Muon; iDecayMode <= kThreeProng1Pi0; ++iDecayMode ) {

      if ( !(iDecayMode == kElectron_Muon  ||
	     iDecayMode == kOneProng0Pi0   ||
	     iDecayMode == kOneProngGt0Pi0 ||
	     iDecayMode == kThreeProng0Pi0) ) continue;
      
      TString decayMode_string = getDecayMode_string(iDecayMode);
      
      double momMin = plotEntries_.front()->momMin_;
      double momMax = plotEntries_.back()->momMax_;
      double momStepSize = 5.;
      unsigned iPad = 0;
      for ( double mom = momMin; mom <= momMax; mom += momStepSize ) {
	if ( (iPad % 6) == 0 ) {
	  canvas->Clear();
	  canvas->Divide(2,3);
	}
	
	TVirtualPad* pad = canvas->cd((iPad % 6) + 1);
	pad->SetLogy();
	
	plotEntryOneMomBin* plotEntry = findPlotEntry(mom);
	if ( plotEntry ) plotEntry->Draw(iDecayMode, option);
	
	if ( (iPad % 6) == 5 ) {
	  canvas->Update();
	  TString outputFileName_i = outputFileName;
	  outputFileName_i = outputFileName_i.ReplaceAll(".", Form("_%s_page%u.", decayMode_string.Data(), (iPad / 6) + 1));
	  canvas->SaveAs(outputFileName_i.Data());
	  iPad = 0;
	} else {
	  ++iPad;
	}
      }
    }
  }

  void Write(TDirectory* dir)
  {
    for ( unsigned iDecayMode = kElectron_Muon; iDecayMode <= kThreeProng1Pi0; ++iDecayMode ) {
      TString decayMode_string = getDecayMode_string(iDecayMode);

      dir->cd();
      TDirectory* subdir = dir->mkdir(decayMode_string);
      subdir->cd();

      for ( std::vector<plotEntryOneMomBin*>::iterator plotEntry = plotEntries_.begin();
	    plotEntry != plotEntries_.end(); ++plotEntry ) {
	(*plotEntry)->getHistogram(iDecayMode)->Write();
      }
    }
  }

  std::vector<plotEntryOneMomBin*> plotEntries_;
};

void Draw(TCanvas* canvas, plotEntryAllMomBins* plots_all, plotEntryAllMomBins* plots_selected1, plotEntryAllMomBins* plots_selected2, 
	  const TString& outputFileName, const char* option)
{
  for ( unsigned iDecayMode = kElectron_Muon; iDecayMode <= kThreeProng1Pi0; ++iDecayMode ) {

    if ( !(iDecayMode == kElectron_Muon  ||
	   iDecayMode == kOneProng0Pi0   ||
	   iDecayMode == kOneProngGt0Pi0 ||
	   iDecayMode == kThreeProng0Pi0) ) continue;

    TString decayMode_string = getDecayMode_string(iDecayMode);

    double momMin = plots_all->plotEntries_.front()->momMin_;
    double momMax = plots_all->plotEntries_.back()->momMax_;
    double momStepSize = 5.;
    unsigned iPad = 0;
    for ( double mom = momMin; mom <= momMax; mom += momStepSize ) {
      if ( (iPad % 6) == 0 ) {
	canvas->Clear();
	canvas->Divide(2,3);
      }

      TVirtualPad* pad = canvas->cd((iPad % 6) + 1);
      pad->SetLogy();

      plotEntryOneMomBin* plotEntry_all       = plots_all->findPlotEntry(mom);
      plotEntryOneMomBin* plotEntry_selected1 = plots_selected1->findPlotEntry(mom);
      plotEntryOneMomBin* plotEntry_selected2 = plots_selected2->findPlotEntry(mom);
      if ( plotEntry_all && plotEntry_selected1 && plotEntry_selected2 ) {
	plotEntry_all->Draw(iDecayMode, option, 1);
	plotEntry_selected1->Draw(iDecayMode, TString(option).Append("same"), 4);
	plotEntry_selected2->Draw(iDecayMode, TString(option).Append("same"), 2);
      }

      if ( (iPad % 6) == 5 ) {
	canvas->Update();
	TString outputFileName_i = outputFileName;
	outputFileName_i = outputFileName_i.ReplaceAll(".", Form("_%s_page%u.", decayMode_string.Data(), (iPad / 6) + 1));
	canvas->SaveAs(outputFileName_i.Data());
	iPad = 0;
      } else {
	++iPad;
      }
    }
  }
}

int main(int argc, const char* argv[])
{
  printTimeStamp("<makeTauDecayKinePlots::main (begin)>");

  TString inputFileNames = "/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*_ggH.root";
  //TString inputFileNames = "/data2/friis/PtBalanceNtupleData_v4/ptBalanceData_mass_200_ggAH.root";

  // CV: chose common binning starting at 15 GeV for all tau decay modes
  TArrayD momBinning_all = getBinningMom("Electron_Muon", "all");
  TArrayD momBinning_selected = getBinningMom("Electron_Muon", "selected");
  TArrayD sepBinning = getBinningSep();
  TArrayD sepTimesMomBinning = getBinningSepTimesMom();

//--- stop RooT from keeping references to all histograms
  TH1::AddDirectory(false); 

  plotEntryAllMomBins* plotsGJDeltaRLab_all = 
    new plotEntryAllMomBins("GJDeltaRLab", "all", momBinning_all, sepBinning);
  plotEntryAllMomBins* plotsGJAngleLab_all = 
    new plotEntryAllMomBins("GJAngleLab", "all", momBinning_all, sepBinning);
  plotEntryAllMomBins* plotsGJDeltaRLab_selected1 = 
    new plotEntryAllMomBins("GJDeltaRLab", "selected1", momBinning_selected, sepBinning);
  plotEntryAllMomBins* plotsGJAngleLab_selected1 = 
    new plotEntryAllMomBins("GJAngleLab", "selected1", momBinning_selected, sepBinning);
  plotEntryAllMomBins* plotsGJDeltaRLab_selected2 = 
    new plotEntryAllMomBins("GJDeltaRLab", "selected2", momBinning_selected, sepBinning);
  plotEntryAllMomBins* plotsGJAngleLab_selected2 = 
    new plotEntryAllMomBins("GJAngleLab", "selected2", momBinning_selected, sepBinning);

  plotEntryAllMomBins* plotsGJDeltaRLabTimesPt_all = 
    new plotEntryAllMomBins("GJDeltaRLabTimesPt", "all", momBinning_all, sepBinning);
  plotEntryAllMomBins* plotsGJAngleLabTimesEnergy_all = 
    new plotEntryAllMomBins("GJAngleLabTimesEnergy", "all", momBinning_all, sepBinning);
  plotEntryAllMomBins* plotsGJDeltaRLabTimesPt_selected1 = 
    new plotEntryAllMomBins("GJDeltaRLabTimesPt", "selected1", momBinning_selected, sepBinning);
  plotEntryAllMomBins* plotsGJAngleLabTimesEnergy_selected1 = 
    new plotEntryAllMomBins("GJAngleLabTimesEnergy", "selected1", momBinning_selected, sepBinning);
  plotEntryAllMomBins* plotsGJDeltaRLabTimesPt_selected2 = 
    new plotEntryAllMomBins("GJDeltaRLabTimesPt", "selected2", momBinning_selected, sepBinning);
  plotEntryAllMomBins* plotsGJAngleLabTimesEnergy_selected2 = 
    new plotEntryAllMomBins("GJAngleLabTimesEnergy", "selected2", momBinning_selected, sepBinning);

  plotEntryAllMomBins* plotsVisInvisDeltaRLab_all = 
    new plotEntryAllMomBins("VisInvisDeltaRLab", "all", momBinning_all, sepBinning);
  plotEntryAllMomBins* plotsVisInvisAngleLab_all = 
    new plotEntryAllMomBins("VisInvisAngleLab", "all", momBinning_all, sepBinning);
  plotEntryAllMomBins* plotsVisInvisDeltaRLab_selected1 = 
    new plotEntryAllMomBins("VisInvisDeltaRLab", "selected1", momBinning_selected, sepBinning);
  plotEntryAllMomBins* plotsVisInvisAngleLab_selected1 = 
    new plotEntryAllMomBins("VisInvisAngleLab", "selected1", momBinning_selected, sepBinning);
  plotEntryAllMomBins* plotsVisInvisDeltaRLab_selected2 = 
    new plotEntryAllMomBins("VisInvisDeltaRLab", "selected2", momBinning_selected, sepBinning);
  plotEntryAllMomBins* plotsVisInvisAngleLab_selected2 = 
    new plotEntryAllMomBins("VisInvisAngleLab", "selected2", momBinning_selected, sepBinning);

  plotEntryAllMomBins* plotsVisInvisDeltaRLabTimesPt_all = 
    new plotEntryAllMomBins("VisInvisDeltaRLabTimesPt", "all", momBinning_all, sepBinning);
  plotEntryAllMomBins* plotsVisInvisAngleLabTimesEnergy_all = 
    new plotEntryAllMomBins("VisInvisAngleLabTimesEnergy", "all", momBinning_all, sepBinning);
  plotEntryAllMomBins* plotsVisInvisDeltaRLabTimesPt_selected1 = 
    new plotEntryAllMomBins("VisInvisDeltaRLabTimesPt", "selected1", momBinning_selected, sepBinning);
  plotEntryAllMomBins* plotsVisInvisAngleLabTimesEnergy_selected1 = 
    new plotEntryAllMomBins("VisInvisAngleLabTimesEnergy", "selected1", momBinning_selected, sepBinning);
  plotEntryAllMomBins* plotsVisInvisDeltaRLabTimesPt_selected2 = 
    new plotEntryAllMomBins("VisInvisDeltaRLabTimesPt", "selected2", momBinning_selected, sepBinning);
  plotEntryAllMomBins* plotsVisInvisAngleLabTimesEnergy_selected2 = 
    new plotEntryAllMomBins("VisInvisAngleLabTimesEnergy", "selected2", momBinning_selected, sepBinning);

  gROOT->SetBatch(true);

  // Load the data
  std::cout << "Loading data from " << inputFileNames <<  std::endl;
  TChain* dataTree = new TChain("makePtBalanceNtuple/ptBalanceNtuple");
  dataTree->Add(inputFileNames);

  /************************************************************************
   ********       Define variabes and data                             ****
   ************************************************************************/

  unsigned numEntries = dataTree->GetEntries();

  std::cout << "Processing " << numEntries << " TTree entries..." << std::endl;

  Double_t leg1Energy, leg1Pt, leg1Eta, leg1VisPt, leg1VisEta, leg1DecayMode;
  Double_t leg1GJDeltaRLab, leg1GJAngleLab, leg1VisInvisAngleLab, leg1VisInvisDeltaRLab;
  dataTree->SetBranchAddress("leg1Energy", &leg1Energy);
  dataTree->SetBranchAddress("leg1Pt", &leg1Pt);
  dataTree->SetBranchAddress("leg1Eta", &leg1Eta);
  dataTree->SetBranchAddress("leg1VisPt", &leg1VisPt);
  dataTree->SetBranchAddress("leg1VisEta", &leg1VisEta);
  dataTree->SetBranchAddress("leg1DecayMode", &leg1DecayMode);
  dataTree->SetBranchAddress("leg1GJDeltaRLab", &leg1GJDeltaRLab);
  dataTree->SetBranchAddress("leg1GJAngleLab", &leg1GJAngleLab);
  dataTree->SetBranchAddress("leg1VisInvisAngleLab", &leg1VisInvisAngleLab);
  dataTree->SetBranchAddress("leg1VisInvisDeltaRLab", &leg1VisInvisDeltaRLab);

  Double_t leg2Energy, leg2Pt, leg2Eta,leg2VisPt, leg2VisEta, leg2DecayMode;
  Double_t leg2GJDeltaRLab, leg2GJAngleLab, leg2VisInvisAngleLab, leg2VisInvisDeltaRLab;
  dataTree->SetBranchAddress("leg2Energy", &leg2Energy);
  dataTree->SetBranchAddress("leg2Pt", &leg2Pt);
  dataTree->SetBranchAddress("leg2Eta", &leg2Eta);
  dataTree->SetBranchAddress("leg2VisPt", &leg2VisPt);
  dataTree->SetBranchAddress("leg2VisEta", &leg2VisEta);
  dataTree->SetBranchAddress("leg2DecayMode", &leg2DecayMode);
  dataTree->SetBranchAddress("leg2GJDeltaRLab", &leg2GJDeltaRLab);
  dataTree->SetBranchAddress("leg2GJAngleLab", &leg2GJAngleLab);
  dataTree->SetBranchAddress("leg2VisInvisAngleLab", &leg2VisInvisAngleLab);
  dataTree->SetBranchAddress("leg2VisInvisDeltaRLab", &leg2VisInvisDeltaRLab);

  for ( unsigned iEntry = 0; iEntry < numEntries; ++iEntry ) {
    dataTree->GetEvent(iEntry); 

    int leg1DecayMode_decoded = decodeDecayMode(leg1DecayMode);
    bool leg1DecayIsHadronic = (leg1DecayMode_decoded >= kOneProng0Pi0 && leg1DecayMode_decoded <= kThreeProng1Pi0);
    int leg2DecayMode_decoded = decodeDecayMode(leg2DecayMode);
    bool leg2DecayIsHadronic = (leg2DecayMode_decoded >= kOneProng0Pi0 && leg2DecayMode_decoded <= kThreeProng1Pi0);
    
    plotsGJDeltaRLab_all->fill(leg1Pt, leg1GJDeltaRLab, leg1DecayMode_decoded);
    plotsGJAngleLab_all->fill(leg1Energy, leg1GJAngleLab, leg1DecayMode_decoded);
    plotsGJDeltaRLabTimesPt_all->fill(leg1Pt, leg1GJDeltaRLab*leg1Pt, leg1DecayMode_decoded);
    plotsGJAngleLabTimesEnergy_all->fill(leg1Energy, leg1GJAngleLab*leg1Energy, leg1DecayMode_decoded);

    plotsGJDeltaRLab_all->fill(leg2Pt, leg2GJDeltaRLab, leg2DecayMode_decoded);
    plotsGJAngleLab_all->fill(leg2Energy, leg2GJAngleLab, leg2DecayMode_decoded);
    plotsGJDeltaRLabTimesPt_all->fill(leg2Pt, leg2GJDeltaRLab*leg2Pt, leg2DecayMode_decoded);
    plotsGJAngleLabTimesEnergy_all->fill(leg2Energy, leg2GJAngleLab*leg2Energy, leg2DecayMode_decoded);

    plotsVisInvisDeltaRLab_all->fill(leg1Pt, leg1VisInvisDeltaRLab, leg1DecayMode_decoded);
    plotsVisInvisAngleLab_all->fill(leg1Energy, leg1VisInvisAngleLab, leg1DecayMode_decoded);
    plotsVisInvisDeltaRLabTimesPt_all->fill(leg1Pt, leg1VisInvisDeltaRLab*leg1Pt, leg1DecayMode_decoded);
    plotsVisInvisAngleLabTimesEnergy_all->fill(leg1Energy, leg1VisInvisAngleLab*leg1Energy, leg1DecayMode_decoded);
    
    plotsVisInvisDeltaRLab_all->fill(leg2Pt, leg2VisInvisDeltaRLab, leg2DecayMode_decoded);
    plotsVisInvisAngleLab_all->fill(leg2Energy, leg2VisInvisAngleLab, leg2DecayMode_decoded);
    plotsVisInvisDeltaRLabTimesPt_all->fill(leg2Pt, leg2VisInvisDeltaRLab*leg2Pt, leg2DecayMode_decoded);
    plotsVisInvisAngleLabTimesEnergy_all->fill(leg2Energy, leg2VisInvisAngleLab*leg2Energy, leg2DecayMode_decoded);

    if ( (leg1DecayMode_decoded == kElectron_Muon && leg1VisPt > 15. && TMath::Abs(leg1VisEta) < 2.1) ||
    	 (leg1DecayIsHadronic                     && leg1VisPt > 20. && TMath::Abs(leg1VisEta) < 2.3) ) {
      plotsGJDeltaRLab_selected1->fill(leg1Pt, leg1GJDeltaRLab, leg1DecayMode_decoded);
      plotsGJAngleLab_selected1->fill(leg1Energy, leg1GJAngleLab, leg1DecayMode_decoded);
      plotsGJDeltaRLabTimesPt_selected1->fill(leg1Pt, leg1GJDeltaRLab*leg1Pt, leg1DecayMode_decoded);
      plotsGJAngleLabTimesEnergy_selected1->fill(leg1Energy, leg1GJAngleLab*leg1Energy, leg1DecayMode_decoded);
    
      plotsVisInvisDeltaRLab_selected1->fill(leg1Pt, leg1VisInvisDeltaRLab, leg1DecayMode_decoded);
      plotsVisInvisAngleLab_selected1->fill(leg1Energy, leg1VisInvisAngleLab, leg1DecayMode_decoded);
      plotsVisInvisDeltaRLabTimesPt_selected1->fill(leg1Pt, leg1VisInvisDeltaRLab*leg1Pt, leg1DecayMode_decoded);
      plotsVisInvisAngleLabTimesEnergy_selected1->fill(leg1Energy, leg1VisInvisAngleLab*leg1Energy, leg1DecayMode_decoded);
    }
    
    if ( (leg2DecayMode_decoded == kElectron_Muon && leg2VisPt > 15. && TMath::Abs(leg2VisEta) < 2.1) ||
    	 (leg2DecayIsHadronic                     && leg2VisPt > 20. && TMath::Abs(leg2VisEta) < 2.3) ) {
      plotsGJDeltaRLab_selected1->fill(leg2Pt, leg2GJDeltaRLab, leg2DecayMode_decoded);
      plotsGJAngleLab_selected1->fill(leg2Energy, leg2GJAngleLab, leg2DecayMode_decoded);
      plotsGJDeltaRLabTimesPt_selected1->fill(leg2Pt, leg2GJDeltaRLab*leg2Pt, leg2DecayMode_decoded);
      plotsGJAngleLabTimesEnergy_selected1->fill(leg2Energy, leg2GJAngleLab*leg2Energy, leg2DecayMode_decoded);

      plotsVisInvisDeltaRLab_selected1->fill(leg2Pt, leg2VisInvisDeltaRLab, leg2DecayMode_decoded);
      plotsVisInvisAngleLab_selected1->fill(leg2Energy, leg2VisInvisAngleLab, leg2DecayMode_decoded);
      plotsVisInvisDeltaRLabTimesPt_selected1->fill(leg2Pt, leg2VisInvisDeltaRLab*leg2Pt, leg2DecayMode_decoded);
      plotsVisInvisAngleLabTimesEnergy_selected1->fill(leg2Energy, leg2VisInvisAngleLab*leg2Energy, leg2DecayMode_decoded);
    }

    if ( ((leg1DecayMode_decoded == kElectron_Muon && leg1VisPt > 15. && TMath::Abs(leg1VisEta) < 2.1) ||
	  (leg1DecayIsHadronic                     && leg1VisPt > 20. && TMath::Abs(leg1VisEta) < 2.3)) &&
	 ((leg2DecayMode_decoded == kElectron_Muon && leg2VisPt > 15. && TMath::Abs(leg2VisEta) < 2.1) ||
	  (leg2DecayIsHadronic                     && leg2VisPt > 20. && TMath::Abs(leg2VisEta) < 2.3)) ) {
      plotsGJDeltaRLab_selected2->fill(leg1Pt, leg1GJDeltaRLab, leg1DecayMode_decoded);
      plotsGJAngleLab_selected2->fill(leg1Energy, leg1GJAngleLab, leg1DecayMode_decoded);
      plotsGJDeltaRLabTimesPt_selected2->fill(leg1Pt, leg1GJDeltaRLab*leg1Pt, leg1DecayMode_decoded);
      plotsGJAngleLabTimesEnergy_selected2->fill(leg1Energy, leg1GJAngleLab*leg1Energy, leg1DecayMode_decoded);
      
      plotsVisInvisDeltaRLab_selected2->fill(leg1Pt, leg1VisInvisDeltaRLab, leg1DecayMode_decoded);
      plotsVisInvisAngleLab_selected2->fill(leg1Energy, leg1VisInvisAngleLab, leg1DecayMode_decoded);
      plotsVisInvisDeltaRLabTimesPt_selected2->fill(leg1Pt, leg1VisInvisDeltaRLab*leg1Pt, leg1DecayMode_decoded);
      plotsVisInvisAngleLabTimesEnergy_selected2->fill(leg1Energy, leg1VisInvisAngleLab*leg1Energy, leg1DecayMode_decoded);

      plotsGJDeltaRLab_selected2->fill(leg2Pt, leg2GJDeltaRLab, leg2DecayMode_decoded);
      plotsGJAngleLab_selected2->fill(leg2Energy, leg2GJAngleLab, leg2DecayMode_decoded);
      plotsGJDeltaRLabTimesPt_selected2->fill(leg2Pt, leg2GJDeltaRLab*leg2Pt, leg2DecayMode_decoded);
      plotsGJAngleLabTimesEnergy_selected2->fill(leg2Energy, leg2GJAngleLab*leg2Energy, leg2DecayMode_decoded);
    
      plotsVisInvisDeltaRLab_selected2->fill(leg2Pt, leg2VisInvisDeltaRLab, leg2DecayMode_decoded);
      plotsVisInvisAngleLab_selected2->fill(leg2Energy, leg2VisInvisAngleLab, leg2DecayMode_decoded);
      plotsVisInvisDeltaRLabTimesPt_selected2->fill(leg2Pt, leg2VisInvisDeltaRLab*leg2Pt, leg2DecayMode_decoded);
      plotsVisInvisAngleLabTimesEnergy_selected2->fill(leg2Energy, leg2VisInvisAngleLab*leg2Energy, leg2DecayMode_decoded);
    }
  }

//--- write histograms to ROOT file
//   (needs to be done before histograms are plotted,
//    as plotting causes histograms to get normalized to unit area)
  TFile* outputFile = new TFile("makeTauDecayKinePlots.root", "RECREATE");

  TDirectory* dirVisInvisDeltaRLab = outputFile->mkdir("VisInvisDeltaRLab");  
  TDirectory* dirVisInvisDeltaRLab_all = dirVisInvisDeltaRLab->mkdir("all");
  plotsVisInvisDeltaRLab_all->Write(dirVisInvisDeltaRLab_all);
  dirVisInvisDeltaRLab->cd();
  TDirectory* dirVisInvisDeltaRLab_selected1 = dirVisInvisDeltaRLab->mkdir("selected1");
  plotsVisInvisDeltaRLab_selected1->Write(dirVisInvisDeltaRLab_selected1);
  outputFile->cd();
  TDirectory* dirVisInvisDeltaRLab_selected2 = dirVisInvisDeltaRLab->mkdir("selected2");
  plotsVisInvisDeltaRLab_selected2->Write(dirVisInvisDeltaRLab_selected2);
  outputFile->cd();

  TDirectory* dirVisInvisAngleLab = outputFile->mkdir("VisInvisAngleLab");
  TDirectory* dirVisInvisAngleLab_all = dirVisInvisAngleLab->mkdir("all");
  plotsVisInvisAngleLab_all->Write(dirVisInvisAngleLab_all);
  dirVisInvisAngleLab->cd();
  TDirectory* dirVisInvisAngleLab_selected1 = dirVisInvisAngleLab->mkdir("selected1");
  plotsVisInvisAngleLab_selected1->Write(dirVisInvisAngleLab_selected1);
  dirVisInvisAngleLab->cd();
  TDirectory* dirVisInvisAngleLab_selected2 = dirVisInvisAngleLab->mkdir("selected2");
  plotsVisInvisAngleLab_selected2->Write(dirVisInvisAngleLab_selected2);

  TDirectory* dirVisInvisDeltaRLabTimesPt = outputFile->mkdir("VisInvisDeltaRLabTimesPt");  
  TDirectory* dirVisInvisDeltaRLabTimesPt_all = dirVisInvisDeltaRLabTimesPt->mkdir("all");
  plotsVisInvisDeltaRLabTimesPt_all->Write(dirVisInvisDeltaRLabTimesPt_all);
  dirVisInvisDeltaRLabTimesPt->cd();
  TDirectory* dirVisInvisDeltaRLabTimesPt_selected1 = dirVisInvisDeltaRLabTimesPt->mkdir("selected1");
  plotsVisInvisDeltaRLabTimesPt_selected1->Write(dirVisInvisDeltaRLabTimesPt_selected1);
  outputFile->cd();
  TDirectory* dirVisInvisDeltaRLabTimesPt_selected2 = dirVisInvisDeltaRLabTimesPt->mkdir("selected2");
  plotsVisInvisDeltaRLabTimesPt_selected2->Write(dirVisInvisDeltaRLabTimesPt_selected2);
  outputFile->cd();

  TDirectory* dirVisInvisAngleLabTimesEnergy = outputFile->mkdir("VisInvisAngleLabTimesEnergy");
  TDirectory* dirVisInvisAngleLabTimesEnergy_all = dirVisInvisAngleLabTimesEnergy->mkdir("all");
  plotsVisInvisAngleLabTimesEnergy_all->Write(dirVisInvisAngleLabTimesEnergy_all);
  dirVisInvisAngleLabTimesEnergy->cd();
  TDirectory* dirVisInvisAngleLabTimesEnergy_selected1 = dirVisInvisAngleLabTimesEnergy->mkdir("selected1");
  plotsVisInvisAngleLabTimesEnergy_selected1->Write(dirVisInvisAngleLabTimesEnergy_selected1);
  dirVisInvisAngleLabTimesEnergy->cd();
  TDirectory* dirVisInvisAngleLabTimesEnergy_selected2 = dirVisInvisAngleLabTimesEnergy->mkdir("selected2");
  plotsVisInvisAngleLabTimesEnergy_selected2->Write(dirVisInvisAngleLabTimesEnergy_selected2);

  delete outputFile;

//--- plot histograms
  TCanvas* canvas = new TCanvas("canvas", "canvas", 1, 1, 600, 900);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
/*
  plotsGJDeltaRLab_all->Draw(canvas, "plots/plotsGJDeltaRLab_all.eps", "e1p");
  plotsGJAngleLab_all->Draw(canvas, "plots/plotsGJAngleLab_all.eps", "e1p");
  plotsGJDeltaRLab_selected1->Draw(canvas, "plots/plotsGJDeltaRLab_selected1.eps", "e1p");
  plotsGJAngleLab_selected1->Draw(canvas, "plots/plotsGJAngleRLab_selected1.eps", "e1p");
  //plotsGJDeltaRLab_selected2->Draw(canvas, "plots/plotsGJDeltaRLab_selected2.eps", "e1p");
  //plotsGJAngleLab_selected2->Draw(canvas, "plots/plotsGJAngleRLab_selected2.eps", "e1p");

  Draw(canvas, plotsGJDeltaRLab_all,
       plotsGJDeltaRLab_selected1, plotsGJDeltaRLab_selected2, 
       "plots/plotsGJDeltaRLab_all_vs_selected.eps", "hist");
  Draw(canvas, plotsGJAngleLab_all, 
       plotsGJAngleLab_selected1, plotsGJDeltaRLab_selected2, 
       "plots/plotsGJAngleLab_all_vs_selected.eps", "hist");

  plotsGJDeltaRLabTimesPt_all->Draw(canvas, "plots/plotsGJDeltaRLabTimesPt_all.eps", "e1p");
  plotsGJAngleLabTimesEnergy_all->Draw(canvas, "plots/plotsGJAngleLabTimesEnergy_all.eps", "e1p");
  plotsGJDeltaRLabTimesPt_selected1->Draw(canvas, "plots/plotsGJDeltaRLabTimesPt_selected1.eps", "e1p");
  plotsGJAngleLabTimesEnergy_selected1->Draw(canvas, "plots/plotsGJAngleRLabTimesEnergy_selected1.eps", "e1p");
  //plotsGJDeltaRLabTimesPt_selected2->Draw(canvas, "plots/plotsGJDeltaRLabTimesPt_selected2.eps", "e1p");
  //plotsGJAngleLabTimesEnergy_selected2->Draw(canvas, "plots/plotsGJAngleRLabTimesEnergy_selected2.eps", "e1p");

  Draw(canvas, plotsGJDeltaRLabTimesPt_all, 
       plotsGJDeltaRLabTimesPt_selected1, plotsGJDeltaRLabTimesPt_selected2, 
       "plots/plotsGJDeltaRLabTimesPt_all_vs_selected.eps", "hist");
  Draw(canvas, plotsGJAngleLabTimesEnergy_all, 
       plotsGJAngleLabTimesEnergy_selected1, plotsGJAngleLabTimesEnergy_selected2,
       "plots/plotsGJAngleLabTimesEnergy_all_vs_selected.eps", "hist");
 */
  plotsVisInvisDeltaRLab_all->Draw(canvas, "plots/plotsVisInvisDeltaRLab_all.eps", "e1p");
  plotsVisInvisAngleLab_all->Draw(canvas, "plots/plotsVisInvisAngleLab_all.eps", "e1p");
  plotsVisInvisDeltaRLab_selected1->Draw(canvas, "plots/plotsVisInvisDeltaRLab_selected1.eps", "e1p");
  plotsVisInvisAngleLab_selected1->Draw(canvas, "plots/plotsVisInvisAngleRLab_selected1.eps", "e1p");
  //plotsVisInvisDeltaRLab_selected2->Draw(canvas, "plots/plotsVisInvisDeltaRLab_selected2.eps", "e1p");
  //plotsVisInvisAngleLab_selected2->Draw(canvas, "plots/plotsVisInvisAngleRLab_selected2.eps", "e1p");

  Draw(canvas, plotsVisInvisDeltaRLab_all, 
       plotsVisInvisDeltaRLab_selected1, plotsVisInvisDeltaRLab_selected2, 
       "plots/plotsVisInvisDeltaRLab_all_vs_selected.eps", "hist");
  Draw(canvas, plotsVisInvisAngleLab_all, 
       plotsVisInvisAngleLab_selected1, plotsVisInvisAngleLab_selected2, 
       "plots/plotsVisInvisAngleLab_all_vs_selected.eps", "hist");

  plotsVisInvisDeltaRLabTimesPt_all->Draw(canvas, "plots/plotsVisInvisDeltaRLabTimesPt_all.eps", "e1p");
  plotsVisInvisAngleLabTimesEnergy_all->Draw(canvas, "plots/plotsVisInvisAngleLabTimesEnergy_all.eps", "e1p");
  plotsVisInvisDeltaRLabTimesPt_selected1->Draw(canvas, "plots/plotsVisInvisDeltaRLabTimesPt_selected1.eps", "e1p");
  plotsVisInvisAngleLabTimesEnergy_selected1->Draw(canvas, "plots/plotsVisInvisAngleRLabTimesEnergy_selected1.eps", "e1p");
  //plotsVisInvisDeltaRLabTimesPt_selected2->Draw(canvas, "plots/plotsVisInvisDeltaRLabTimesPt_selected2.eps", "e1p");
  //plotsVisInvisAngleLabTimesEnergy_selected2->Draw(canvas, "plots/plotsVisInvisAngleRLabTimesEnergy_selected2.eps", "e1p");

  Draw(canvas, plotsVisInvisDeltaRLabTimesPt_all, 
       plotsVisInvisDeltaRLabTimesPt_selected1, plotsVisInvisDeltaRLabTimesPt_selected2, 
       "plots/plotsVisInvisDeltaRLabTimesPt_all_vs_selected.eps", "hist");
  Draw(canvas, plotsVisInvisAngleLabTimesEnergy_all, 
       plotsVisInvisAngleLabTimesEnergy_selected1, plotsVisInvisAngleLabTimesEnergy_selected2, 
       "plots/plotsVisInvisAngleLabTimesEnergy_all_vs_selected.eps", "hist");

  delete dataTree;

  delete plotsGJDeltaRLab_all;
  delete plotsGJAngleLab_all;
  delete plotsGJDeltaRLab_selected1;
  delete plotsGJAngleLab_selected1;
  delete plotsGJDeltaRLab_selected2;
  delete plotsGJAngleLab_selected2;

  delete plotsGJDeltaRLabTimesPt_all;
  delete plotsGJAngleLabTimesEnergy_all;
  delete plotsGJDeltaRLabTimesPt_selected1;
  delete plotsGJAngleLabTimesEnergy_selected1;
  delete plotsGJDeltaRLabTimesPt_selected2;
  delete plotsGJAngleLabTimesEnergy_selected2;

  delete plotsVisInvisDeltaRLab_all;
  delete plotsVisInvisAngleLab_all;
  delete plotsVisInvisDeltaRLab_selected1;
  delete plotsVisInvisAngleLab_selected1;
  delete plotsVisInvisDeltaRLab_selected2;
  delete plotsVisInvisAngleLab_selected2;

  delete plotsVisInvisDeltaRLabTimesPt_all;
  delete plotsVisInvisAngleLabTimesEnergy_all;
  delete plotsVisInvisDeltaRLabTimesPt_selected1;
  delete plotsVisInvisAngleLabTimesEnergy_selected1;
  delete plotsVisInvisDeltaRLabTimesPt_selected2;
  delete plotsVisInvisAngleLabTimesEnergy_selected2;

  printTimeStamp("<makeTauDecayKinePlots::main (end)>");
}
