#ifndef TauAnalysis_GenSimTools_tauDecayKineAuxFunctions_h
#define TauAnalysis_GenSimTools_tauDecayKineAuxFunctions_h

#include <TString.h>
#include <TArrayD.h>
#include <TMath.h>

enum { kElectron_Muon, 
       kOneProng0Pi0, kOneProng1Pi0, kOneProng2Pi0, kOneProngGt0Pi0, 
       kThreeProng0Pi0, kThreeProng1Pi0 };

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

double encodeDecayMode(int decayMode)
{
  double decayMode_encoded = -1000;
  if      ( decayMode == kElectron_Muon  ) decayMode_encoded = -1.;
  else if ( decayMode == kOneProng0Pi0   ) decayMode_encoded =  0.;
  else if ( decayMode == kOneProng1Pi0   ) decayMode_encoded =  1.;
  else if ( decayMode == kOneProng2Pi0   ) decayMode_encoded =  2.;
  else if ( decayMode == kThreeProng0Pi0 ) decayMode_encoded = 10.;
  else if ( decayMode == kThreeProng1Pi0 ) decayMode_encoded = 11.;
  return decayMode_encoded;
}

//
//-------------------------------------------------------------------------------
//

TArrayD getBinningMom()
{
  TArrayD binning;
  
  double xMin =  15.0;
  double xMax = 350.0;
  double dx   =   1.0;

  unsigned numBins = TMath::Nint((xMax - xMin)/dx);

  binning.Set(numBins + 1);

  for ( unsigned iBin = 0; iBin <= numBins; ++iBin ) {
    double x = xMin + iBin*dx;
    binning[iBin] = x;
  }
  
  return binning;
}

TArrayD getBinningSep()
{
  unsigned numBins = 400;
  double xMin = 0.0;
  double xMax = 1.0;

  TArrayD binning(numBins + 1);
  
  double dx   = (xMax - xMin)/numBins;

  for ( unsigned iBin = 0; iBin <= numBins; ++iBin ) {
    double x = xMin + iBin*dx;
    binning[iBin] = x;
  }
  
  return binning;
}

TArrayD getBinningSepTimesMom()
{
  unsigned numBins = 120;
  double xMin =  0.0;
  double xMax = 12.0;

  TArrayD binning(numBins + 1);
  
  double dx   = (xMax - xMin)/numBins;
  
  for ( unsigned iBin = 0; iBin <= numBins; ++iBin ) {
    double x = xMin + iBin*dx;
    binning[iBin] = x;
  }
  
  return binning;
}

#endif
