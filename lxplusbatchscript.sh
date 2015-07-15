#!/bin/bash
cd /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/
#cmsenv
eval `scramv1 runtime -sh`

cd /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/

nJobs=0
i=0
while [ $i -le $nJobs ];
do 
#   let "start=i*2"
#   let "end=(i+1)*2"
	let "start=0"
	let "end=0"
  export FIRST=$start  
  export LAST=$end 
  echo "First = $FIRST and last file = $LAST"   
  bsub -R "pool>5000" -M 3000000 -q 1nd -J merge_job_${i} < /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/submit.sh
  let "i++"
done

echo "submit all jobs!"

#echo "Copying output files to " $destination