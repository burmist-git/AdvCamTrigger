#!/bin/sh

#simHomeDir="/home/users/b/burmistr/pyeventio_example/"
simHomeDir="./"
n_data_chunks=20

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d    : default - simulation with terzinag4"
    echo " [1]       : particle type (g,gd,e,p)"
    echo " [2]       : njobs"
    echo " [3]       : jobType (sbatch | screen | live)"
    echo " [0] -NGB  : NGB"
    echo " [1]       : njobs"
    echo " [2]       : jobType (sbatch | screen | live)"
    echo " [0] -info : print info"
    echo " [0] -kill : kill all jobs"
    echo " [0] -h    : print help"
}

counter=0

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	if [ $# -eq 4 ]; then
	    particletype=$2
	    njobs=$3
	    typeJob=$4
	    #
	    Ebin=0
	    Thbin=0
	    rbin=0
	    ((njobs=njobs-1))
	    #
	    if [ "$particletype" = "g" ]; then
                nDir=10
                n_multiply=25
            elif [ "$particletype" = "gd" ]; then
                nDir=20
                n_multiply=75
            elif [ "$particletype" = "e" ]; then
                nDir=20
                n_multiply=10
            elif [ "$particletype" = "p" ]; then
                nDir=50
                n_multiply=30
            fi
	    #
	    for i in $(seq 0 $njobs)
	    do
		((counter=counter+1))
		jobID=`printf "%04d" $i`
		echo "jobID = $jobID"
		data_chunk_ID=$(echo "$i/$nDir" | bc)
		#echo "$data_chunk_ID"
		#$simHomeDir/run_trgA_job.sh -d $particletype $Ebin $Thbin $rbin $jobID $data_chunk_ID
		if [ "$typeJob" = "sbatch" ]; then
		    sbatch $simHomeDir/run_trgA_job.sh -d $particletype $Ebin $Thbin $rbin $jobID $data_chunk_ID $typeJob
		elif [ "$typeJob" = "screen" ]; then
		    $simHomeDir/run_trgA_job.sh -d $particletype $Ebin $Thbin $rbin $jobID $data_chunk_ID $typeJob
		elif [ "$typeJob" = "live" ]; then
		    $simHomeDir/run_trgA_job.sh -d $particletype $Ebin $Thbin $rbin $jobID $data_chunk_ID $typeJob
		else
		    printHelp
		fi
	    done
	else
	    printHelp
	fi
    elif [ "$1" = "-NGB" ]; then
	if [ $# -eq 3 ]; then
	    njobs=$2
	    typeJob=$3
	    ((njobs=njobs-1))
	    for i in $(seq 0 $njobs)
	    do
		jobID=`printf "%04d" $i`
		echo "jobID = $jobID"
		#
		if [ "$typeJob" = "sbatch" ]; then
		    sbatch $simHomeDir/run_trgA_job.sh -NGB $jobID $typeJob
		elif [ "$typeJob" = "screen" ]; then
		    $simHomeDir/run_trgA_job.sh -NGB $jobID $typeJob
		elif [ "$typeJob" = "live" ]; then
		    $simHomeDir/run_trgA_job.sh -NGB $jobID $typeJob
		else
		    printHelp
		fi
	    done	
	else
	    printHelp
	fi	
    elif [ "$1" = "-info" ]; then
	squeue | head -n 1
	squeue | grep burmistr
    elif [ "$1" = "-kill" ]; then
	scancel --user=burmistr --state=pending
	scancel --user=burmistr --state=CG
	scancel --user=burmistr --state=R
    elif [ "$1" = "-h" ]; then
	printHelp
    else
        printHelp
    fi
fi

#espeak "I have done"
