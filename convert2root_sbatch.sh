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
	#
	#
	i_start=0
	i_stop=19
	for i in `seq $i_start $i_stop`
        do
            folderID=`printf %04d $i`
	    sbatch /home/users/b/burmistr/pyeventio_example/convert2root_job.sh -d gd $folderID $i
	    #/home/users/b/burmistr/pyeventio_example/convert2root_job.sh -d gd $folderID $i
        done
	#
	#
	i_start=0
	i_stop=49
	for i in `seq $i_start $i_stop`
        do
            folderID=`printf %04d $i`
	    sbatch /home/users/b/burmistr/pyeventio_example/convert2root_job.sh -d p $folderID $i
	    #/home/users/b/burmistr/pyeventio_example/convert2root_job.sh -d p $folderID $i
        done
	#
	#
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
