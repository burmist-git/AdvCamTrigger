#!/bin/bash

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d  : build apptainer (singularity) no modules need to be loaded"
    echo " [0] -h  : print help"
}


if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	#
	singularity --version
	#
	rm -rf pyeventio_example_singularity.sif
	singularity build pyeventio_example_singularity.sif pyeventio_example_singularity.def
	#
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
