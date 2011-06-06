#!/bin/csh

#echo "starting Angle vs. Energy fits..."
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' Electron_Muon leg1 angle all >&! fitTauDecayKine_Electron_Muon_angle_all.out &
nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' Electron_Muon leg1 angle selected >&! fitTauDecayKine_Electron_Muon_angle_selected.out &
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' OneProng0Pi0 leg2 angle all >&! fitTauDecayKine_OneProng0Pi0_angle_all.out &
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' OneProng0Pi0 leg2 angle selected >&! fitTauDecayKine_OneProng0Pi0_angle_selected.out &
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' OneProngGt0Pi0 leg2 angle all >&! fitTauDecayKine_OneProngGt0Pi0_angle_all.out &
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' OneProngGt0Pi0 leg2 angle selected >&! fitTauDecayKine_OneProngGt0Pi0_angle_selected.out &
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' ThreeProng0Pi0 leg2 angle all >&! fitTauDecayKine_ThreeProng0Pi0_angle_all.out &
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' ThreeProng0Pi0 leg2 angle selected >&! fitTauDecayKine_ThreeProng0Pi0_angle_selected.out &

#echo "starting dR vs. Pt fits..."
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' Electron_Muon leg1 dR all >&! fitTauDecayKine_Electron_Muon_dR_all.out &
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' Electron_Muon leg1 dR selected >&! fitTauDecayKine_Electron_Muon_dR_selected.out &
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' OneProng0Pi0 leg2 dR all >&! fitTauDecayKine_OneProng0Pi0_dR_all.out &
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' OneProng0Pi0 leg2 dR selected >&! fitTauDecayKine_OneProng0Pi0_dR_selected.out &
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' OneProngGt0Pi0 leg2 dR all >&! fitTauDecayKine_OneProngGt0Pi0_dR_all.out &
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' OneProngGt0Pi0 leg2 dR selected >&! fitTauDecayKine_OneProngGt0Pi0_dR_selected.out &
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' ThreeProng0Pi0 leg2 dR all >&! fitTauDecayKine_ThreeProng0Pi0_dR_all.out &
#nice ../../../../bin/slc5_amd64_gcc434/fitTauDecayKine '/data2/friis/PtBalanceNtupleData_v5/ptBalanceData_*ggH.root' ThreeProng0Pi0 leg2 dR selected >&! fitTauDecayKine_ThreeProng0Pi0_dR_selected.out &

