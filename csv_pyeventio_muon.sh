#!/bin/bash

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -m : muon"
    echo " [0] -h : print help"
}

if [ $# -eq 0 ]; then
    printHelp
else
    if [ "$1" = "-m" ]; then
	#
        datafilein="../scratch/simtel_data/muon/data/run1_muon.simtel.gz" 
	headerout="../scratch/simtel_data/muon/csv/run1_muon.header.csv"
        pe_info_out="../scratch/simtel_data/muon/csv/run1_muon.pe_info.csv"
	#
	rm -rf $headerout
        rm -rf $pe_info_out
	#
	python csv_pyeventio.py $datafilein $headerout $pe_info_out
	#
    elif [ "$1" = "-h" ]; then
	printHelp
    else
        printHelp
    fi
fi
