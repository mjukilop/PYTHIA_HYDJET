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
	fi

	if [ $2 -eq 15 ]
	then
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
	else
		echo "Arguments were $1, $2"
		echo "Argument one is the number of Jobs you wish to run per file, two is either a specific pT file (e.g. 15,30,80,270) or "all" to run over all files"
	fi
fi

echo "submit all jobs!"

#echo "Copying output files to " $destination
