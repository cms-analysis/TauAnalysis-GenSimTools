#!/bin/bash

export SCRAM_ARCH=slc5_ia32_gcc434
export VO_CMS_SW_DIR=/afs/cern.ch/cms/
source $VO_CMS_SW_DIR/cmsset_default.sh

cd $LS_SUBCWD
eval `scram runtime -sh`
cd -
# Append the working directory to arguments

cmsRun $LS_SUBCWD/$@
export result=$?
#mv -f *.root $LS_SUBCWD
rfcp *.root /castor/cern.ch/cms/store/user/friis/NSVfitNtuple/
exit $result
