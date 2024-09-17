#!/bin/bash

# Wed Sep 11 04:30:20 PM CEST 2024
# Autor: Leonid Burmistrov

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] --test              : (test)"
    echo " [0] --inc_reduce_stack  : increase reduce stack memory"
    echo " [0] --set_unlimited_mem : set unlimited memory"
    echo " [0] -h                  : print help"
}
    
if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "--test" ]; then
	make clean; make;
	header_file="header.csv"
	pe_info_file_LST1="pe_info_LST1.csv"
	pe_info_file_LST2="pe_info_LST2.csv"
	pe_info_file_LST3="pe_info_LST3.csv"
	pe_info_file_LST4="pe_info_LST4.csv"
 	outputRootFile="corsika_proton_run1.root"
	#
	./convert2rootstereo 0 $header_file $pe_info_file_LST1 $pe_info_file_LST2 $pe_info_file_LST3 $pe_info_file_LST4 $outputRootFile
    elif [ "$1" = "--inc_reduce_stack" ]; then	
	echo " ulimit -s           : current limit"
	echo " ulimit -s unlimited : set stack unlimited"
	echo " sudo swapoff -a; sudo swapon -a : clean swap"
    elif [ "$1" = "--set_unlimited_mem" ]; then	
	ulimit -s unlimited
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
