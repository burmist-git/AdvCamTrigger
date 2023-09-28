#!/bin/sh

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d    : default"
    echo " [0] -info : print info"
    echo " [0] -kill : kill all jobs"
    echo " [0] -h    : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	conda activate pyeventio
	for i in $(ls /srv/beegfs/scratch/users/b/burmistr/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_diffuse_nsb_1x/pkl/); do
	    jobID=$i
	    sbatch /home/users/b/burmistr/pyeventio_example/pkl_tocsv_job.sh -d gd $jobID
	done
	for i in $(ls /srv/beegfs/scratch/users/b/burmistr/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/pkl/); do
	    jobID=$i
	    sbatch /home/users/b/burmistr/pyeventio_example/pkl_tocsv_job.sh -d p $jobID
	done
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
