#!/usr/bin/env python

import os
import socket
import subprocess
import time

inputFileNames = "/castor/cern.ch/user/v/veelken/CMSSW_4_1_x/ntuples/tauDecayKine/v5/ptBalanceData_*_ggH.root"

decayModesToRun = [
    "Electron_Muon",
    "OneProng0Pi0",
    "OneProngGt0Pi0",
    "ThreeProng0Pi0"	
]

legsToRun = [
    "leg1",
    "leg2"	
]

parametrizationsToRun = [
    "dR",
    "angle"	
]

selectionsToRun = [
    "all",
    "selected"	
]	

submit = "yes"
#submit = "no"

# use '1nw' queue for "full" fit, '1nh' queue for prefit only
#queue = "1nw"
queue = "1nh"

#resourceRequest = None
resourceRequest = "rusage[mem=3200, swp=3200]"

workingDirectory = os.getcwd()
architecture = os.environ["SCRAM_ARCH"]
hostname = "ucdavis.cern.ch"
plotsDirectory = "/data1/veelken/CMSSW_4_1_x/plots/GenSimTools"

# iterate over all combinations of decay modes and seletions
for decayMode in decayModesToRun:
    for leg in legsToRun:
        for parametrization in parametrizationsToRun:
            for selection in selectionsToRun:

                inputFilePath = inputFileNames[:inputFileNames.rfind('/')]
    	        print("inputFilePath = %s" % inputFilePath)
    	        inputFileName_pattern = inputFileNames[inputFileNames.rfind('/') + 1:]
	        print("inputFileName_pattern = %s" % inputFileName_pattern)

	        # build shell script
	        submissionDirectory = os.path.join(workingDirectory, "lxbatch")
	        print("submissionDirectory = %s" % submissionDirectory)
	        if not os.path.exists(submissionDirectory):
	            os.makedirs(submissionDirectory)
        
                print("plotsDirectory = %s" % plotsDirectory)

	        executable = "%s/bin/%s/fitTauDecayKine" % (workingDirectory[:workingDirectory.find("/src/")], architecture)
	        print("executable = %s" % executable)
	        parameter = "'%s' %s %s %s %s" % (inputFileName_pattern, decayMode, leg, parametrization, selection)
	        print("parameter = %s" % parameter)

	        nslsCommand = "nsls %s" % inputFilePath	
	        inputFileName_pattern_items = inputFileName_pattern.split('*')
	        for inputFileName_pattern_item in inputFileName_pattern_items:
	            nslsCommand += " | grep %s" % inputFileName_pattern_item
	        print("nslsCommand = %s" % nslsCommand)

	        script = """#!/bin/csh
limit vmem unlim
cd %(submissionDirectory)s
setenv SCRAM_ARCH %(architecture)s
eval `scram runtime -csh`
cd -
#set fileNames=(`%(nslsCommand)s`)
#foreach fileName (${fileNames})
#    echo "copying fileName %(inputFilePath)s/${fileName} --> ./${fileName}"
#    rfcp %(inputFilePath)s/${fileName} ./${fileName}
#end
rfcp /castor/cern.ch/user/v/veelken/CMSSW_4_1_x/plots/tauDecayKine/makeTauDecayKinePlots.root .
mkdir plots
%(executable)s %(parameter)s
scp ./fitTauDecayKinePlots*.root %(workingDirectory)s
scp ./fitTauDecayKinePlots*.txt %(workingDirectory)s
scp ./plots/* %(hostname)s:%(plotsDirectory)s
scp ./mcTauDecayKine*.root %(workingDirectory)s
""" % {
    'submissionDirectory' : submissionDirectory,
    'architecture' : architecture,
    'nslsCommand' : nslsCommand,
    'inputFilePath' : inputFilePath,
    'executable' : executable,
    'parameter' :  parameter,
    'workingDirectory' : workingDirectory,
    'hostname' : hostname,
    'plotsDirectory' : plotsDirectory
}
	        scriptFileName = os.path.join(submissionDirectory, "runFitTauDecayKine_%s_%s_%s_%s.csh" \
                                % (decayMode, leg, parametrization, selection))
	        scriptFile = open(scriptFileName, "w")
	        scriptFile.write(script)
	        scriptFile.close()    

                # make shell script executable
	        os.chmod(scriptFileName, 0744)
    
	        logFileDirectory = os.path.join(workingDirectory, "lxbatch_log")
	        if not os.path.exists(logFileDirectory):
	            os.makedirs(logFileDirectory)
    
                # submit job to CERN batch system
	        if submit == "yes":
	            logFileName = os.path.join(logFileDirectory, "fitTauDecayKine_%s_%s_%s_%s.log" \
                                 % (decayMode, leg, parametrization, selection))
	            print("logFileName = " + logFileName)
	            jobName = "fitTauDecayKine_%s_%s_%s_%s" % (decayMode, leg, parametrization, selection)
                    bsubCommand = 'bsub -q ' + queue + ' -J ' + jobName + ' -L /bin/csh -eo ' + logFileName + ' -oo ' + logFileName
	            if resourceRequest != None:
	                bsubCommand += ' -R \"' + resourceRequest + '\" '			
	            bsubCommand += ' < ' + scriptFileName
	            print("bsubCommand = '%s'" % bsubCommand)
	            subprocess.call(bsubCommand, shell = True)

	            # wait for a few seconds, in order not to generate too many castor requests in too short a time
	            # (maximum number of permissible requests = 900 in 180 seconds; jobs will abort in case limit gets exceeded)
  	            time.sleep(1)
    
