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
        pe_info_out_LST1="pe_info_LST1.csv"
	pe_info_out_LST2="pe_info_LST2.csv"
	pe_info_out_LST3="pe_info_LST3.csv"
	pe_info_out_LST4="pe_info_LST4.csv"
	#
	python csv_pyeventio_stereo.py $datafilein $headerout $pe_info_out_LST1 $pe_info_out_LST2 $pe_info_out_LST3 $pe_info_out_LST4
	#
    elif [ "$1" = "-h" ]; then
	printHelp
    else
        printHelp
    fi
fi
