#!/bin/csh

../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v4/ptBalanceData_*ggAH.root' Electron_Muon leg1 dR all >&! fitTauDecayKine_Electron_Muon.out
#../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v4/ptBalanceData_*ggAH.root' OneProng0Pi0 leg2 dR all >&! fitTauDecayKine_OneProng0Pi0.out
