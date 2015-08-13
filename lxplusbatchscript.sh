#!/bin/bash
cd /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/
#cmsenv
eval `scramv1 runtime -sh`

cd /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/

if [ $# -eq 2 ]
then
	nJobs=$1
	if [ $2 = "all" ]
	then
		for j in 1 2 3 4 5 6 7 8 9 10
		do
			i=0
			while [ $i -lt $nJobs ]
			do 
				let "job=i"
				let "file=j"
				export FILE=$file
				export JOB=$job 
				export JOBS=$nJobs
				echo "File = $file"
				echo "#Jobs = $nJobs"   
				bsub -R "pool>5000" -g /skeeton -M 3000000 -q 1nd -J merge_job_${i} < /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/submit.sh
				let "i++"
			done
		done
	elif [ $2 -eq 15 ]
	then
		i=0
		while [ $i -lt $nJobs ]
		do 
			let "job=i"
			let "file=1"
			export FILE=$file
			export JOB=$job 
			export JOBS=$nJobs
			echo "File = $file"
			echo "#Jobs = $nJobs" 
			bsub -R "pool>5000" -g /skeeton -M 3000000 -q 1nd -J merge_job_${i} < /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/submit.sh
			let "i++"
		done
	elif [ $2 -eq 30 ]
	then
		i=0
		while [ $i -lt $nJobs ]
		do 
			let "job=i"
			let "file=2"
			export FILE=$file
			export JOB=$job 
			export JOBS=$nJobs
			echo "File = $file"
			echo "#Jobs = $nJobs"  
			bsub -R "pool>5000" -g /skeeton -M 3000000 -q 1nd -J merge_job_${i} < /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/submit.sh
			let "i++"
		done
	elif [ $2 -eq 50 ]
	then
		i=0
		while [ $i -lt $nJobs ]
		do 
			let "job=i"
			let "file=3"
			export FILE=$file
			export JOB=$job 
			export JOBS=$nJobs
			echo "File = $file"
			echo "#Jobs = $nJobs"  
			bsub -R "pool>5000" -g /skeeton -M 3000000 -q 1nd -J merge_job_${i} < /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/submit.sh
			let "i++"
		done
	elif [ $2 -eq 80 ]
	then
		i=0
		while [ $i -lt $nJobs ]
		do 
			let "job=i"
			let "file=4"
			export FILE=$file
			export JOB=$job 
			export JOBS=$nJobs
			echo "File = $file"
			echo "#Jobs = $nJobs"  
			bsub -R "pool>5000" -g /skeeton -M 3000000 -q 1nd -J merge_job_${i} < /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/submit.sh
			let "i++"
		done
	elif [ $2 -eq 120 ]
	then
		i=0
		while [ $i -lt $nJobs ]
		do 
			let "job=i"
			let "file=5"
			export FILE=$file
			export JOB=$job 
			export JOBS=$nJobs
			echo "File = $file"
			echo "#Jobs = $nJobs"  
			bsub -R "pool>5000" -g /skeeton -M 3000000 -q 1nd -J merge_job_${i} < /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/submit.sh
			let "i++"
		done
	elif [ $2 -eq 170 ]
	then
		i=0
		while [ $i -lt $nJobs ]
		do 
			let "job=i"
			let "file=6"
			export FILE=$file
			export JOB=$job 
			export JOBS=$nJobs
			echo "File = $file"
			echo "#Jobs = $nJobs"  
			bsub -R "pool>5000" -g /skeeton -M 3000000 -q 1nd -J merge_job_${i} < /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/submit.sh
			let "i++"
		done
	elif [ $2 -eq 220 ]
	then
		i=0
		while [ $i -lt $nJobs ]
		do 
			let "job=i"
			let "file=7"
			export FILE=$file
			export JOB=$job 
			export JOBS=$nJobs
			echo "File = $file"
			echo "#Jobs = $nJobs"  
			bsub -R "pool>5000" -g /skeeton -M 3000000 -q 1nd -J merge_job_${i} < /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/submit.sh
			let "i++"
		done
	elif [ $2 -eq 280 ]
	then
		i=0
		while [ $i -lt $nJobs ]
		do 
			let "job=i"
			let "file=8"
			export FILE=$file
			export JOB=$job 
			export JOBS=$nJobs
			echo "File = $file"
			echo "#Jobs = $nJobs"  
			bsub -R "pool>5000" -g /skeeton -M 3000000 -q 1nd -J merge_job_${i} < /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/submit.sh
			let "i++"
		done
	elif [ $2 -eq 370 ]
	then
		i=0
		while [ $i -lt $nJobs ]
		do 
			let "job=i"
			let "file=10"
			export FILE=$file
			export JOB=$job 
			export JOBS=$nJobs
			echo "File = $file"
			echo "#Jobs = $nJobs"  
			bsub -R "pool>5000" -g /skeeton -M 3000000 -q 1nd -J merge_job_${i} < /afs/cern.ch/user/s/skeeton/myanalyser/CMSSW_5_3_20/src/root/submit.sh
			let "i++"
		done
	fi
fi

echo "Arguments were $1, $2"
echo "Expected '"lxplusbatchscript nJobsPerFile file"', where file=all,15,30,50,etc"

echo "submit all jobs!"

#echo "Copying output files to " $destination
