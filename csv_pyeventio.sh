#!/bin/bash

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d : default/test"
    echo " [0] -h : print help"
}

if [ $# -eq 0 ]; then
    printHelp
else
    if [ "$1" = "-d" ]; then
	#
        datafilein="/scratch/snx3000/lburmist/simtel_data/proton/data/corsika_run1.simtel.gz" 
        headerout="header.csv"
        pe_info_out="pe_info.csv"
	#
	python csv_pyeventio.py $datafilein $headerout $pe_info_out
	#
    elif [ "$1" = "-h" ]; then
	printHelp
    else
        printHelp
    fi
fi
