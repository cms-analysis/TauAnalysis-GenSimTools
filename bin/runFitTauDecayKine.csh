#!/bin/csh

#../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' Electron_Muon leg1 angle all >&! fitTauDecayKine_Electron_Muon.out
#../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' Electron_Muon leg1 angle selected >&! fitTauDecayKine_Electron_Muon.out
../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' OneProng0Pi0 leg2 angle selected >&! fitTauDecayKine_OneProng0Pi0.out
#../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' ThreeProng0Pi0 leg2 angle all >&! fitTauDecayKine_ThreeProng0Pi0.out
