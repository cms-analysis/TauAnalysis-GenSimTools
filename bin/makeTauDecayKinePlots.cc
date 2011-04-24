
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

enum { kElectron_Muon, 
       kOneProng0Pi0, kOneProng1Pi0, kOneProng2Pi0, kOneProngGt0Pi0, 
       kThreeProng0Pi0, kThreeProng1Pi0 };

TH1* bookHistogram(const TString& decayMode, const TString& variable, const TString& label, double ptMin, double ptMax)
{
  TString histogramName = Form("%s_%sPt%2.0fto%2.0f", decayMode.Data(), label.Data(), ptMin, ptMax);
  TString histogramTitle = Form("%s %s %s: Pt%2.0f..%2.0f", variable.Data(), decayMode.Data(), label.Data(), ptMin, ptMax);
  TH1* histogram = new TH1F(histogramName.Data(), histogramTitle.Data(), 50, 0., 0.5);
  return histogram;
}

TString getDecayMode_string(int decayMode)
{
  TString decayMode_string = "undefined";
  if      ( decayMode == kElectron_Muon  ) decayMode_string = "Electron_Muon";
  else if ( decayMode == kOneProng0Pi0   ) decayMode_string = "OneProng0Pi0";
  else if ( decayMode == kOneProng1Pi0   ) decayMode_string = "OneProng1Pi0"; 
  else if ( decayMode == kOneProng2Pi0   ) decayMode_string = "OneProng2Pi0"; 
  else if ( decayMode == kOneProngGt0Pi0 ) decayMode_string = "OneProngGt0Pi0";
  else if ( decayMode == kThreeProng0Pi0 ) decayMode_string = "ThreeProng0Pi0";
  else if ( decayMode == kThreeProng1Pi0 ) decayMode_string = "ThreeProng1Pi0";
  return decayMode_string;
}

struct plotEntryOnePtBin
{
  plotEntryOnePtBin(const TString& variable, const TString& label, double ptMin, double ptMax)
    : ptMin_(ptMin),
      ptMax_(ptMax)
  {
    deltaElectron_Muon_ = bookHistogram(getDecayMode_string(kElectron_Muon), variable, label, ptMin, ptMax);
    deltaOneProng0Pi0_ = bookHistogram(getDecayMode_string(kOneProng0Pi0), variable, label, ptMin, ptMax);
    deltaOneProng1Pi0_ = bookHistogram(getDecayMode_string(kOneProng1Pi0), variable, label, ptMin, ptMax);
    deltaOneProng2Pi0_ = bookHistogram(getDecayMode_string(kOneProng2Pi0), variable, label, ptMin, ptMax);
    deltaOneProngGt0Pi0_ = bookHistogram(getDecayMode_string(kOneProngGt0Pi0), variable, label, ptMin, ptMax);
    deltaThreeProng0Pi0_ = bookHistogram(getDecayMode_string(kThreeProng0Pi0), variable, label, ptMin, ptMax);
    deltaThreeProng1Pi0_ = bookHistogram(getDecayMode_string(kThreeProng1Pi0), variable, label, ptMin, ptMax);
  }

  ~plotEntryOnePtBin()
  {
    delete deltaElectron_Muon_;
    delete deltaOneProng0Pi0_;
    delete deltaOneProng1Pi0_;
    delete deltaOneProng2Pi0_;
    delete deltaOneProngGt0Pi0_;
    delete deltaThreeProng0Pi0_;
    delete deltaThreeProng1Pi0_;
  }

  void fill(double pt, double delta, int decayMode)
  {
    if ( pt > ptMin_ && pt < ptMax_ ) {
      TH1* histogram = getHistogram(decayMode);
      if ( histogram ) histogram->Fill(delta);
      if ( decayMode == kOneProng1Pi0 || decayMode == kOneProng2Pi0 ) deltaOneProngGt0Pi0_->Fill(delta);
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
      histogram->SetMarkerStyle(20);
      histogram->SetMarkerSize(1.0);
      histogram->SetLineColor(color);
      histogram->SetLineWidth(2);
      histogram->Draw(option);
    }
  }

  TH1* getHistogram(int decayMode)
  {
    if      ( decayMode == kElectron_Muon  ) return deltaElectron_Muon_;
    else if ( decayMode == kOneProng0Pi0   ) return deltaOneProng0Pi0_;
    else if ( decayMode == kOneProng1Pi0   ) return deltaOneProng1Pi0_;
    else if ( decayMode == kOneProng2Pi0   ) return deltaOneProng2Pi0_;    
    else if ( decayMode == kOneProngGt0Pi0 ) return deltaOneProngGt0Pi0_;
    else if ( decayMode == kThreeProng0Pi0 ) return deltaThreeProng0Pi0_;
    else if ( decayMode == kThreeProng1Pi0 ) return deltaThreeProng1Pi0_;
    else return 0;
  }

  double ptMin_;
  double ptMax_;

  TH1* deltaElectron_Muon_;
  TH1* deltaOneProng0Pi0_;
  TH1* deltaOneProng1Pi0_;
  TH1* deltaOneProng2Pi0_;
  TH1* deltaOneProngGt0Pi0_;
  TH1* deltaThreeProng0Pi0_;
  TH1* deltaThreeProng1Pi0_;
};

struct plotEntryAllPtBins
{
  plotEntryAllPtBins(const TString& variable, const TString& label, unsigned numPtBins, const double* ptBinning)
  {
    for ( unsigned iPtBin = 0; iPtBin < numPtBins; ++iPtBin ) {
      double ptMin = ptBinning[iPtBin];
      double ptMax = ptBinning[iPtBin + 1];
      plotEntryOnePtBin* newPlotEntry = new plotEntryOnePtBin(variable, label, ptMin, ptMax);
      plotEntries_.push_back(newPlotEntry);
    }
  }
  
  ~plotEntryAllPtBins()
  {
    for ( std::vector<plotEntryOnePtBin*>::iterator it = plotEntries_.begin();
	  it != plotEntries_.end(); ++it ) {
      delete (*it);
    }
  }

  void fill(double pt, double delta, int genDecayMode)
  {
    for ( std::vector<plotEntryOnePtBin*>::iterator plotEntry = plotEntries_.begin();
	  plotEntry != plotEntries_.end(); ++plotEntry ) {
      (*plotEntry)->fill(pt, delta, genDecayMode);
    }
  }

  void Draw(TCanvas* canvas, const TString& outputFileName, const char* option)
  {
    for ( unsigned iDecayMode = kElectron_Muon; iDecayMode <= kThreeProng1Pi0; ++iDecayMode ) {
      TString decayMode_string = getDecayMode_string(iDecayMode);
      
      unsigned numPtBins = plotEntries_.size();
      for ( unsigned iPtBin = 0; iPtBin < numPtBins; ++iPtBin ) {
	if ( (iPtBin % 6) == 0 ) {
	  canvas->Clear();
	  canvas->Divide(2,3);
	}
	
	TVirtualPad* pad = canvas->cd((iPtBin % 6) + 1);
	pad->SetLogy();
	plotEntries_[iPtBin]->Draw(iDecayMode, option);
	
	if ( (iPtBin % 6) == 5 ) {
	  canvas->Update();
	  TString outputFileName_i = outputFileName;
	  outputFileName_i = outputFileName_i.ReplaceAll(".", Form("_%s_page%u.", decayMode_string.Data(), (iPtBin / 6) + 1));
	  canvas->SaveAs(outputFileName_i.Data());
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

      for ( std::vector<plotEntryOnePtBin*>::iterator plotEntry = plotEntries_.begin();
	    plotEntry != plotEntries_.end(); ++plotEntry ) {
	(*plotEntry)->getHistogram(iDecayMode)->Write();
      }
    }
  }

  std::vector<plotEntryOnePtBin*> plotEntries_;
};

bool equals(double value, double target)
{
  return (TMath::Abs(value - target) < 1.e-2);
}

int decodeDecayMode(double decayMode)
{
  int decayMode_decoded = -1;
  if      ( equals(decayMode, -1) ) decayMode_decoded = kElectron_Muon;
  else if ( equals(decayMode,  0) ) decayMode_decoded = kOneProng0Pi0;
  else if ( equals(decayMode,  1) ) decayMode_decoded = kOneProng1Pi0;
  else if ( equals(decayMode,  2) ) decayMode_decoded = kOneProng2Pi0;
  else if ( equals(decayMode, 10) ) decayMode_decoded = kThreeProng0Pi0;
  else if ( equals(decayMode, 11) ) decayMode_decoded = kThreeProng1Pi0;
  return decayMode_decoded;
}

void Draw(TCanvas* canvas, plotEntryAllPtBins* plots_all, plotEntryAllPtBins* plots_selected1, plotEntryAllPtBins* plots_selected2, 
	  const TString& outputFileName, const char* option)
{
  for ( unsigned iDecayMode = kElectron_Muon; iDecayMode <= kThreeProng1Pi0; ++iDecayMode ) {
    TString decayMode_string = getDecayMode_string(iDecayMode);

    unsigned numPtBins = plots_all->plotEntries_.size();
    for ( unsigned iPtBin = 0; iPtBin < numPtBins; ++iPtBin ) {
      if ( (iPtBin % 6) == 0 ) {
	canvas->Clear();
	canvas->Divide(2,3);
      }

      TVirtualPad* pad = canvas->cd((iPtBin % 6) + 1);
      pad->SetLogy();
      plots_all->plotEntries_[iPtBin]->Draw(iDecayMode, option, 1);
      plots_selected1->plotEntries_[iPtBin]->Draw(iDecayMode, TString(option).Append("same"), 4);
      plots_selected2->plotEntries_[iPtBin]->Draw(iDecayMode, TString(option).Append("same"), 2);

      if ( (iPtBin % 6) == 5 ) {
	canvas->Update();
	TString outputFileName_i = outputFileName;
	outputFileName_i = outputFileName_i.ReplaceAll(".", Form("_%s_page%u.", decayMode_string.Data(), (iPtBin / 6) + 1));
	canvas->SaveAs(outputFileName_i.Data());
      }
    }
  }
}

int main(int argc, const char* argv[])
{
  TString inputFileNames = "/data2/friis/PtBalanceNtupleData_v2/ptBalanceData_*.root";

  //unsigned numPtBins = 12;
  //const double ptBinning[] = { 15., 20., 25., 30., 35., 40., 50., 60., 80., 100., 120., 160., 200. };
  unsigned numPtBins = 36;
  const double ptBinning[] = { 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75.,
			       80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135., 140., 
			       145., 150., 155., 160., 165., 170., 175., 180., 185., 190., 195. };

//--- stop RooT from keeping references to all histograms
  TH1::AddDirectory(false); 

  plotEntryAllPtBins* plotsGJDeltaRLab_all = new plotEntryAllPtBins("GJDeltaRLab", "all", numPtBins, ptBinning);
  plotEntryAllPtBins* plotsGJAngleLab_all = new plotEntryAllPtBins("GJAngleLab", "all", numPtBins, ptBinning);
  plotEntryAllPtBins* plotsGJDeltaRLab_selected1 = new plotEntryAllPtBins("GJDeltaRLab", "selected1", numPtBins, ptBinning);
  plotEntryAllPtBins* plotsGJAngleLab_selected1 = new plotEntryAllPtBins("GJAngleLab", "selected1", numPtBins, ptBinning);
  plotEntryAllPtBins* plotsGJDeltaRLab_selected2 = new plotEntryAllPtBins("GJDeltaRLab", "selected2", numPtBins, ptBinning);
  plotEntryAllPtBins* plotsGJAngleLab_selected2 = new plotEntryAllPtBins("GJAngleLab", "selected2", numPtBins, ptBinning);

  plotEntryAllPtBins* plotsVisInvisDeltaRLab_all = new plotEntryAllPtBins("GJDeltaRLab", "all", numPtBins, ptBinning);
  plotEntryAllPtBins* plotsVisInvisAngleLab_all = new plotEntryAllPtBins("GJAngleLab", "all", numPtBins, ptBinning);
  plotEntryAllPtBins* plotsVisInvisDeltaRLab_selected1 = new plotEntryAllPtBins("GJDeltaRLab", "selected1", numPtBins, ptBinning);
  plotEntryAllPtBins* plotsVisInvisAngleLab_selected1 = new plotEntryAllPtBins("GJAngleLab", "selected1", numPtBins, ptBinning);
  plotEntryAllPtBins* plotsVisInvisDeltaRLab_selected2 = new plotEntryAllPtBins("GJDeltaRLab", "selected2", numPtBins, ptBinning);
  plotEntryAllPtBins* plotsVisInvisAngleLab_selected2 = new plotEntryAllPtBins("GJAngleLab", "selected2", numPtBins, ptBinning);  

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
    
    plotsGJDeltaRLab_all->fill(leg2Pt, leg2GJDeltaRLab, leg2DecayMode_decoded);
    plotsGJAngleLab_all->fill(leg2Energy, leg2GJAngleLab, leg2DecayMode_decoded);

    plotsVisInvisDeltaRLab_all->fill(leg1Pt, leg1VisInvisDeltaRLab, leg1DecayMode_decoded);
    plotsVisInvisAngleLab_all->fill(leg1Energy, leg1VisInvisAngleLab, leg1DecayMode_decoded);
    
    plotsVisInvisDeltaRLab_all->fill(leg2Pt, leg2VisInvisDeltaRLab, leg2DecayMode_decoded);
    plotsVisInvisAngleLab_all->fill(leg2Energy, leg2VisInvisAngleLab, leg2DecayMode_decoded);

    if ( (leg1DecayMode_decoded == kElectron_Muon && leg1VisPt > 15. && TMath::Abs(leg1VisEta) < 2.1) ||
    	 (leg1DecayIsHadronic                     && leg1VisPt > 20. && TMath::Abs(leg1VisEta) < 2.3) ) {
      plotsGJDeltaRLab_selected1->fill(leg1Pt, leg1GJDeltaRLab, leg1DecayMode_decoded);
      plotsGJAngleLab_selected1->fill(leg1Energy, leg1GJAngleLab, leg1DecayMode_decoded);
    
      plotsVisInvisDeltaRLab_selected1->fill(leg1Pt, leg1VisInvisDeltaRLab, leg1DecayMode_decoded);
      plotsVisInvisAngleLab_selected1->fill(leg1Energy, leg1VisInvisAngleLab, leg1DecayMode_decoded);
    }
    
    if ( (leg2DecayMode_decoded == kElectron_Muon && leg2VisPt > 15. && TMath::Abs(leg2VisEta) < 2.1) ||
    	 (leg2DecayIsHadronic                     && leg2VisPt > 20. && TMath::Abs(leg2VisEta) < 2.3) ) {
      plotsGJDeltaRLab_selected1->fill(leg2Pt, leg2GJDeltaRLab, leg2DecayMode_decoded);
      plotsGJAngleLab_selected1->fill(leg2Energy, leg2GJAngleLab, leg2DecayMode_decoded);
    
      plotsVisInvisDeltaRLab_selected1->fill(leg2Pt, leg2VisInvisDeltaRLab, leg2DecayMode_decoded);
      plotsVisInvisAngleLab_selected1->fill(leg2Energy, leg2VisInvisAngleLab, leg2DecayMode_decoded);
    }

    if ( ((leg1DecayMode_decoded == kElectron_Muon && leg1VisPt > 15. && TMath::Abs(leg1VisEta) < 2.1) ||
	  (leg1DecayIsHadronic                     && leg1VisPt > 20. && TMath::Abs(leg1VisEta) < 2.3)) &&
	 ((leg2DecayMode_decoded == kElectron_Muon && leg2VisPt > 15. && TMath::Abs(leg2VisEta) < 2.1) ||
	  (leg2DecayIsHadronic                     && leg2VisPt > 20. && TMath::Abs(leg2VisEta) < 2.3)) ) {
      plotsGJDeltaRLab_selected2->fill(leg1Pt, leg1GJDeltaRLab, leg1DecayMode_decoded);
      plotsGJAngleLab_selected2->fill(leg1Energy, leg1GJAngleLab, leg1DecayMode_decoded);
      
      plotsVisInvisDeltaRLab_selected2->fill(leg1Pt, leg1VisInvisDeltaRLab, leg1DecayMode_decoded);
      plotsVisInvisAngleLab_selected2->fill(leg1Energy, leg1VisInvisAngleLab, leg1DecayMode_decoded);

      plotsGJDeltaRLab_selected2->fill(leg2Pt, leg2GJDeltaRLab, leg2DecayMode_decoded);
      plotsGJAngleLab_selected2->fill(leg2Energy, leg2GJAngleLab, leg2DecayMode_decoded);
    
      plotsVisInvisDeltaRLab_selected2->fill(leg2Pt, leg2VisInvisDeltaRLab, leg2DecayMode_decoded);
      plotsVisInvisAngleLab_selected2->fill(leg2Energy, leg2VisInvisAngleLab, leg2DecayMode_decoded);
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
  TDirectory* dirVisInvisDeltaRLab_selected = dirVisInvisDeltaRLab->mkdir("selected2");
  plotsVisInvisDeltaRLab_selected2->Write(dirVisInvisDeltaRLab_selected);
  outputFile->cd();
  TDirectory* dirVisInvisAngleLab = outputFile->mkdir("VisInvisAngleLab");
  TDirectory* dirVisInvisAngleLab_all = dirVisInvisAngleLab->mkdir("all");
  plotsVisInvisAngleLab_all->Write(dirVisInvisAngleLab_all);
  dirVisInvisAngleLab->cd();
  TDirectory* dirVisInvisAngleLab_selected = dirVisInvisAngleLab->mkdir("selected2");
  plotsVisInvisAngleLab_selected2->Write(dirVisInvisAngleLab_selected);
  delete outputFile;

//--- plot histograms
  TCanvas* canvas = new TCanvas("canvas", "canvas", 1, 1, 600, 900);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  plotsGJDeltaRLab_all->Draw(canvas, "plots/plotsGJDeltaRLab_all.eps", "e1p");
  plotsGJAngleLab_all->Draw(canvas, "plots/plotsGJAngleLab_all.eps", "e1p");
  plotsGJDeltaRLab_selected1->Draw(canvas, "plots/plotsGJDeltaRLab_selected1.eps", "e1p");
  plotsGJAngleLab_selected1->Draw(canvas, "plots/plotsGJAngleRLab_selected1.eps", "e1p");
  plotsGJDeltaRLab_selected2->Draw(canvas, "plots/plotsGJDeltaRLab_selected2.eps", "e1p");
  plotsGJAngleLab_selected2->Draw(canvas, "plots/plotsGJAngleRLab_selected2.eps", "e1p");

  Draw(canvas, plotsGJDeltaRLab_all, plotsGJDeltaRLab_selected1, plotsGJDeltaRLab_selected2, 
       "plots/plotsGJDeltaRLab_all_vs_selected.eps", "hist");
  Draw(canvas, plotsGJAngleLab_all, plotsGJAngleLab_selected1, plotsGJDeltaRLab_selected2, 
       "plots/plotsGJAngleLab_all_vs_selected.eps", "hist");

  plotsVisInvisDeltaRLab_all->Draw(canvas, "plots/plotsVisInvisDeltaRLab_all.eps", "e1p");
  plotsVisInvisAngleLab_all->Draw(canvas, "plots/plotsVisInvisAngleLab_all.eps", "e1p");
  plotsVisInvisDeltaRLab_selected1->Draw(canvas, "plots/plotsVisInvisDeltaRLab_selected1.eps", "e1p");
  plotsVisInvisAngleLab_selected1->Draw(canvas, "plots/plotsVisInvisAngleRLab_selected1.eps", "e1p");
  plotsVisInvisDeltaRLab_selected2->Draw(canvas, "plots/plotsVisInvisDeltaRLab_selected2.eps", "e1p");
  plotsVisInvisAngleLab_selected2->Draw(canvas, "plots/plotsVisInvisAngleRLab_selected2.eps", "e1p");

  Draw(canvas, plotsVisInvisDeltaRLab_all, plotsVisInvisDeltaRLab_selected1, plotsVisInvisDeltaRLab_selected2, 
       "plots/plotsVisInvisDeltaRLab_all_vs_selected.eps", "hist");
  Draw(canvas, plotsVisInvisAngleLab_all, plotsVisInvisAngleLab_selected1, plotsVisInvisAngleLab_selected2, 
       "plots/plotsVisInvisAngleLab_all_vs_selected.eps", "hist");

  delete dataTree;

  delete plotsGJDeltaRLab_all;
  delete plotsGJAngleLab_all;
  delete plotsGJDeltaRLab_selected1;
  delete plotsGJAngleLab_selected1;
  delete plotsGJDeltaRLab_selected2;
  delete plotsGJAngleLab_selected2;

  delete plotsVisInvisDeltaRLab_all;
  delete plotsVisInvisAngleLab_all;
  delete plotsVisInvisDeltaRLab_selected1;
  delete plotsVisInvisAngleLab_selected1;
  delete plotsVisInvisDeltaRLab_selected2;
  delete plotsVisInvisAngleLab_selected2;
}
