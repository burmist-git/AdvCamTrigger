#!/bin/bash

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -m                  : muon"
    echo " [0] --inc_reduce_stack  : increase reduce stack memory"
    echo " [0] --set_unlimited_mem : set unlimited memory"
    echo " [0] -h                  : print help"
}
    
if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-m" ]; then
	header_file="../scratch/simtel_data/muon/csv/run1_muon.header.csv"
	pe_info_file="../scratch/simtel_data/muon/csv/run1_muon.pe_info.csv"
	outputRootFile="../scratch/simtel_data/muon/root/run1_muon.root"
	#
	rm -rf $outputRootFile
	#
	./convert2root 1 $header_file $pe_info_file $outputRootFile
    elif [ "$1" = "--inc_reduce_stack" ]; then	
	echo " ulimit -s                       : current limit"
	echo " ulimit -s unlimited             : set stack unlimited"
	echo " sudo swapoff -a; sudo swapon -a : clean swap"
    elif [ "$1" = "--set_unlimited_mem" ]; then	
	ulimit -s unlimited
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
