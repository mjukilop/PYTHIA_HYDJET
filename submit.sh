#!/bin/bash

cd /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/
#cmsenv
eval `scramv1 runtime -sh`


#Added by Ian
export X509_USER_PROXY=~/x509_user_proxy/proxy
voms-proxy-init --noregen
#</Ian>

cd /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/

echo "root -l -b -q macro.c"
echo "I am in $PWD"   

root -b -q macro.c <<EOF
#.x macro.c;
#.q
EOF


echo "Done all jobs!"

destination="/afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/"

echo "Copying output files to " $destination

mv myhisto.root $destination 