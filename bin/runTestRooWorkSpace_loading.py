#!/usr/bin/env python

import os

executable = "../../../../bin/slc5_amd64_gcc434/testRooWorkSpace_loading"

#inputFilePath = "/afs/cern.ch/user/v/veelken/scratch0/CMSSW_4_1_3/src/TauAnalysis/GenSimTools/test/mcTauDecayKine_2011May17"
inputFilePath = "/afs/cern.ch/user/v/veelken/scratch0/CMSSW_4_1_3/src/TauAnalysis/GenSimTools/bin"

wsList = [
    # tau --> e/mu nu nu Decays, dR vs. Pt
    [ "mcTauDecayKine_Electron_Muon_leg1_dR_all_ws_prefit.root", "ws_prefit", 
      "pdf_Electron_Muon_AllMom_leg1VisInvisDeltaRLab_leg1_dR_all", "leg1Pt", "leg1VisInvisDeltaRLabTimesPt" ],
    #[ "mcTauDecayKine_Electron_Muon_leg1_dR_all_ws_fit.root", "ws_fit", 
    #  "pdf_Electron_Muon_AllMom_leg1VisInvisDeltaRLab_leg1_dR_all", "leg1Pt", "leg1VisInvisDeltaRLabTimesPt" ],
    [ "mcTauDecayKine_Electron_Muon_leg1_dR_selected_ws_prefit.root", "ws_prefit", 
      "pdf_Electron_Muon_AllMom_leg1VisInvisDeltaRLab_leg1_dR_selected1", "leg1Pt", "leg1VisInvisDeltaRLabTimesPt" ],
    #[ "mcTauDecayKine_Electron_Muon_leg1_dR_selected_ws_fit.root", "ws_fit", 
    #  "pdf_Electron_Muon_AllMom_leg1VisInvisDeltaRLab_leg1_dR_selected1", "leg1Pt", "leg1VisInvisDeltaRLabTimesPt" ],

    # tau --> e/mu nu nu Decays, 3d-angle vs. Energy
    [ "mcTauDecayKine_Electron_Muon_leg1_angle_all_ws_prefit.root", "ws_prefit", 
      "pdf_Electron_Muon_AllMom_leg1VisInvisAngleLab_leg1_angle_all", "leg1Energy", "leg1VisInvisAngleLabTimesEnergy" ],
    #[ "mcTauDecayKine_Electron_Muon_leg1_angle_all_ws_fit.root", "ws_fit", 
    #  "pdf_Electron_Muon_AllMom_leg1VisInvisAngleLab_leg1_angle_all", "leg1Energy", "leg1VisInvisAngleLabTimesEnergy" ],
    [ "mcTauDecayKine_Electron_Muon_leg1_angle_selected_ws_prefit.root", "ws_prefit", 
      "pdf_Electron_Muon_AllMom_leg1VisInvisAngleLab_leg1_angle_selected1", "leg1Energy", "leg1VisInvisAngleLabTimesEnergy" ],
    #[ "mcTauDecayKine_Electron_Muon_leg1_angle_selected_ws_fit.root", "ws_fit", 
    #  "pdf_Electron_Muon_AllMom_leg1VisInvisAngleLab_leg1_angle_selected1", "leg1Energy", "leg1VisInvisAngleLabTimesEnergy" ],
     
    # tau --> pi- nu Decays, dR vs. Pt
    [ "mcTauDecayKine_OneProng0Pi0_leg2_dR_all_ws_prefit.root", "ws_prefit", 
      "pdf_OneProng0Pi0_AllMom_leg2VisInvisDeltaRLab_leg2_dR_all", "leg2Pt", "leg2VisInvisDeltaRLabTimesPt" ],
    #[ "mcTauDecayKine_OneProng0Pi0_leg2_dR_all_ws_fit.root", "ws_fit", 
    #  "pdf_OneProng0Pi0_AllMom_leg2VisInvisDeltaRLab_leg2_dR_all", "leg2Pt", "leg2VisInvisDeltaRLabTimesPt" ],
    [ "mcTauDecayKine_OneProng0Pi0_leg2_dR_selected_ws_prefit.root", "ws_prefit", 
      "pdf_OneProng0Pi0_AllMom_leg2VisInvisDeltaRLab_leg2_dR_selected1", "leg2Pt", "leg2VisInvisDeltaRLabTimesPt" ],
    #[ "mcTauDecayKine_OneProng0Pi0_leg2_dR_selected_ws_fit.root", "ws_fit", 
    #  "pdf_OneProng0Pi0_AllMom_leg2VisInvisDeltaRLab_leg2_dR_selected1", "leg2Pt", "leg2VisInvisDeltaRLabTimesPt" ],

    # tau --> pi- nu Decays, 3d-angle vs. Energy
    [ "mcTauDecayKine_OneProng0Pi0_leg2_angle_all_ws_prefit.root", "ws_prefit", 
      "pdf_OneProng0Pi0_AllMom_leg2VisInvisAngleLab_leg2_angle_all", "leg2Energy", "leg2VisInvisAngleLabTimesEnergy" ],
    #[ "mcTauDecayKine_OneProng0Pi0_leg2_angle_all_ws_fit.root", "ws_fit", 
    #  "pdf_OneProng0Pi0_AllMom_leg2VisInvisAngleLab_leg2_angle_all", "leg2Energy", "leg2VisInvisAngleLabTimesEnergy" ],
    [ "mcTauDecayKine_OneProng0Pi0_leg2_angle_selected_ws_prefit.root", "ws_prefit", 
      "pdf_OneProng0Pi0_AllMom_leg2VisInvisAngleLab_leg2_angle_selected1", "leg2Energy", "leg2VisInvisAngleLabTimesEnergy" ],
    #[ "mcTauDecayKine_OneProng0Pi0_leg2_angle_selected_ws_fit.root", "ws_fit", 
    #  "pdf_OneProng0Pi0_AllMom_leg2VisInvisAngleLab_leg2_angle_selected1", "leg2Energy", "leg2VisInvisAngleLabTimesEnergy" ],

    # tau --> pi- 1(2) pi0 nu Decays, dR vs. Pt
    [ "mcTauDecayKine_OneProngGt0Pi0_leg2_dR_all_ws_prefit.root", "ws_prefit", 
      "pdf_OneProngGt0Pi0_AllMom_leg2VisInvisDeltaRLab_leg2_dR_all", "leg2Pt", "leg2VisInvisDeltaRLabTimesPt" ],
    #[ "mcTauDecayKine_OneProngGt0Pi0_leg2_dR_all_ws_fit.root", "ws_fit", 
    #  "pdf_OneProngGt0Pi0_AllMom_leg2VisInvisDeltaRLab_leg2_dR_all", "leg2Pt", "leg2VisInvisDeltaRLabTimesPt" ],
    [ "mcTauDecayKine_OneProngGt0Pi0_leg2_dR_selected_ws_prefit.root", "ws_prefit", 
      "pdf_OneProngGt0Pi0_AllMom_leg2VisInvisDeltaRLab_leg2_dR_selected1", "leg2Pt", "leg2VisInvisDeltaRLabTimesPt" ],
    #[ "mcTauDecayKine_OneProngGt0Pi0_leg2_dR_selected_ws_fit.root", "ws_fit", 
    #  "pdf_OneProngGt0Pi0_AllMom_leg2VisInvisDeltaRLab_leg2_dR_selected1", "leg2Pt", "leg2VisInvisDeltaRLabTimesPt" ],

    # tau --> pi- 1(2) pi0 nu Decays, 3d-angle vs. Energy
    [ "mcTauDecayKine_OneProngGt0Pi0_leg2_angle_all_ws_prefit.root", "ws_prefit", 
      "pdf_OneProngGt0Pi0_AllMom_leg2VisInvisAngleLab_leg2_angle_all", "leg2Energy", "leg2VisInvisAngleLabTimesEnergy" ],
    #[ "mcTauDecayKine_OneProngGt0Pi0_leg2_angle_all_ws_fit.root", "ws_fit", 
    #  "pdf_OneProngGt0Pi0_AllMom_leg2VisInvisAngleLab_leg2_angle_all", "leg2Energy", "leg2VisInvisAngleLabTimesEnergy" ],
    [ "mcTauDecayKine_OneProngGt0Pi0_leg2_angle_selected_ws_prefit.root", "ws_prefit", 
      "pdf_OneProngGt0Pi0_AllMom_leg2VisInvisAngleLab_leg2_angle_selected1", "leg2Energy", "leg2VisInvisAngleLabTimesEnergy" ],
    #[ "mcTauDecayKine_OneProngGt0Pi0_leg2_angle_selected_ws_fit.root", "ws_fit", 
    #  "pdf_OneProngGt0Pi0_AllMom_leg2VisInvisAngleLab_leg2_angle_selected1", "leg2Energy", "leg2VisInvisAngleLabTimesEnergy" ],

    # tau --> pi- pi+ pi- nu Decays, dR vs. Pt
    [ "mcTauDecayKine_ThreeProng0Pi0_leg2_dR_all_ws_prefit.root", "ws_prefit", 
      "pdf_ThreeProng0Pi0_AllMom_leg2VisInvisDeltaRLab_leg2_dR_all", "leg2Pt", "leg2VisInvisDeltaRLabTimesPt" ],
    #[ "mcTauDecayKine_ThreeProng0Pi0_leg2_dR_all_ws_fit.root", "ws_fit", 
    #  "pdf_ThreeProng0Pi0_AllMom_leg2VisInvisDeltaRLab_leg2_dR_all", "leg2Pt", "leg2VisInvisDeltaRLabTimesPt" ],
    [ "mcTauDecayKine_ThreeProng0Pi0_leg2_dR_selected_ws_prefit.root", "ws_prefit", 
      "pdf_ThreeProng0Pi0_AllMom_leg2VisInvisDeltaRLab_leg2_dR_selected1", "leg2Pt", "leg2VisInvisDeltaRLabTimesPt" ],
    #[ "mcTauDecayKine_ThreeProng0Pi0_leg2_dR_selected_ws_fit.root", "ws_fit", 
    #  "pdf_ThreeProng0Pi0_AllMom_leg2VisInvisDeltaRLab_leg2_dR_selected1", "leg2Pt", "leg2VisInvisDeltaRLabTimesPt" ],

    # tau --> pi- pi+ pi- nu Decays, 3d-angle vs. Energy
    [ "mcTauDecayKine_ThreeProng0Pi0_leg2_angle_all_ws_prefit.root", "ws_prefit", 
      "pdf_ThreeProng0Pi0_AllMom_leg2VisInvisAngleLab_leg2_angle_all", "leg2Energy", "leg2VisInvisAngleLabTimesEnergy" ],
    #[ "mcTauDecayKine_ThreeProng0Pi0_leg2_angle_all_ws_fit.root", "ws_fit", 
    #  "pdf_ThreeProng0Pi0_AllMom_leg2VisInvisAngleLab_leg2_angle_all", "leg2Energy", "leg2VisInvisAngleLabTimesEnergy" ],
    [ "mcTauDecayKine_ThreeProng0Pi0_leg2_angle_selected_ws_prefit.root", "ws_prefit", 
      "pdf_ThreeProng0Pi0_AllMom_leg2VisInvisAngleLab_leg2_angle_selected1", "leg2Energy", "leg2VisInvisAngleLabTimesEnergy" ],
    #[ "mcTauDecayKine_ThreeProng0Pi0_leg2_angle_selected_ws_fit.root", "ws_fit", 
    #  "pdf_ThreeProng0Pi0_AllMom_leg2VisInvisAngleLab_leg2_angle_selected1", "leg2Energy", "leg2VisInvisAngleLabTimesEnergy" ],
]

shFileName = "runTestRooWorkSpace_loading.csh"

shFileContent = ""
shFileContent += "#!/bin/csh\n"
shFileContent += "\n"

for ws in wsList:

    if len(ws) != 5:
	raise ValueError("RooWorkSpace definition has invalid format !!")

    inputFileName = ws[0]
    wsName = ws[1]
    pdfName = ws[2]
    momName = ws[3]
    sepTimesMomName = ws[4]

    logFileName = inputFileName
    logFileName = logFileName.replace(".root", ".log")

    shFileContent += "%s %s %s %s %s %s >&! %s\n" \
      % (executable, os.path.join(inputFilePath, inputFileName), wsName, pdfName, momName, sepTimesMomName, logFileName)

shFile = open(shFileName, "w")
shFile.write(shFileContent)
shFile.close()

# make shell script executable
os.chmod(shFileName, 0744)

print "\n"
print "Run ./%s to start jobs" % shFileName
