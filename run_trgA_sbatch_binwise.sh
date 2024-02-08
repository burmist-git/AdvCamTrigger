#!/bin/sh

simHomeDir="/home/users/b/burmistr/pyeventio_example/"
#simHomeDir="./"

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d    : default - simulation with terzinag4"
    echo " [1]       : particle type (g,gd,e,p)"
    echo " [0] -NGB  : NGB"
    echo " [0]       : njobs"
    echo " [0] -info : print info"
    echo " [0] -kill : kill all jobs"
    echo " [0] -h    : print help"
}

Ebin_arr=(
    5
    6
    7
    8
    9
    10
    11
    12
    13
    14
    15
    16
    17
    18
    19
)

Thbin_arr=(
    0
    1
    2
    3
    4
    5
    6
    7
    8
    9
)

rbin_arr=(
    0
    1
    2
    3
    4
    5
    6
    7
    8
    9    
)

counter=0

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	if [ $# -eq 2 ]; then
	    particletype=$2
            for Ebin in ${Ebin_arr[@]} ; do
		for Thbin in ${Thbin_arr[@]} ; do
		    for rbin in ${rbin_arr[@]} ; do
			#echo "Ebin  = $Ebin"
			#echo "Thbin = $Thbin"
			#echo "rbin  = $rbin"
			jobID=`printf "%04d" $counter`
			((counter=counter+1))
			echo "jobID = $jobID"
			data_chunk_ID="-999"
			#$simHomeDir/run_trgA_job.sh -d $particletype $Ebin $Thbin $rbin $jobID
			sbatch $simHomeDir/run_trgA_job.sh -d $particletype $Ebin $Thbin $rbin $jobID $data_chunk_ID
		    done
		done
	    done	
	else
	    printHelp
	fi
    elif [ "$1" = "-NGB" ]; then
	if [ $# -eq 2 ]; then
	    njobs=$2
	    ((njobs=njobs-1))
	    for i in $(seq 0 $njobs)
	    do
		jobID=`printf "%04d" $i`
		echo "jobID = $jobID"
		sbatch $simHomeDir/run_trgA_job.sh -NGB $jobID
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
