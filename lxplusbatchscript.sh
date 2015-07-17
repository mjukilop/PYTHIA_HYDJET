#!/bin/bash
cd /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/
#cmsenv
eval `scramv1 runtime -sh`

cd /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/

nJobs=2
i=0
while [ $i -lt $nJobs ];
do 
#   let "start=i*2"
#   let "end=(i+1)*2"
	let "start=0"
	let "end=0"
	let "job=i"
  export FIRST=$start  
  export LAST=$end
  export JOB=$job 
  export JOBS=$nJobs
  echo "Job = $JOB"   
  bsub -R "pool>5000" -M 3000000 -q 8nm -J merge_job_${i} < /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/submit.sh
  let "i++"
done

echo "submit all jobs!"

#echo "Copying output files to " $destination
