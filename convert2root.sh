#!/bin/bash

# Thu 11 May 16:09:40 CEST 2023
# Autor: Leonid Burmistrov



function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] --convert_gm       : convert gammas"
    echo " [0] --convert_gd       : convert gammas diffuse"
    echo " [0] --convert_el       : convert electrons"
    echo " [0] --convert_pr       : convert protons"
    echo " [0] --inc_reduce_stack : increase reduce stack memory "
    echo " [0] -h                 : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "--convert_gm" ]; then
	make clean; make;
	particle="gamma"
	#
	header_file="../compressed_data/no_nsb_cut/$particle/corsika_run307.header.csv"
	pe_info_file="../compressed_data/no_nsb_cut/$particle/corsika_run307.pe_info.csv"
	wf_file_list="./wf_file_list_no_nsb_cut.list"
	outputRootFile="../compressed_data/no_nsb_cut/$particle/corsika_run307.compressed.root"
	#
	#header_file="../compressed_data/$particle/corsika_run307.header.csv"
	#pe_info_file="../compressed_data/$particle/corsika_run307.pe_info.csv"
	#wf_file_list="./wf_file_list.list"
	#outputRootFile="../compressed_data/$particle/corsika_run307.compressed.root"
	#
	./convert2root 0 $header_file $pe_info_file $wf_file_list $outputRootFile
    elif [ "$1" = "--convert_gd" ]; then
	make clean; make;
	particle="gamma_diffuse"
	header_file="../compressed_data/no_nsb_cut/$particle/corsika_run307.header.csv"
	pe_info_file="../compressed_data/no_nsb_cut/$particle/corsika_run307.pe_info.csv"
	outputRootFile="../compressed_data/no_nsb_cut/$particle/corsika_run307.compressed.root"
	#./convert2root 0 $header_file $pe_info_file $outputRootFile
    elif [ "$1" = "--convert_el" ]; then
	make clean; make;
	particle="electron"
	header_file="../compressed_data/no_nsb_cut/$particle/corsika_run307.header.csv"
	pe_info_file="../compressed_data/no_nsb_cut/$particle/corsika_run307.pe_info.csv"
	outputRootFile="../compressed_data/no_nsb_cut/$particle/corsika_run307.compressed.root"
	#./convert2root 0 $header_file $pe_info_file $outputRootFile
    elif [ "$1" = "--convert_pr" ]; then
	make clean; make;
	particle="proton"
	header_file="../compressed_data/no_nsb_cut/$particle/corsika_run307.header.csv"
	pe_info_file="../compressed_data/no_nsb_cut/$particle/corsika_run307.pe_info.csv"
	outputRootFile="../compressed_data/no_nsb_cut/$particle/corsika_run307.compressed.root"
	#./convert2root 0 $header_file $pe_info_file $outputRootFile
    elif [ "$1" = "--inc_reduce_stack" ]; then	
	echo " ulimit -s           : current limit"
	echo " ulimit -s unlimited : set stack unlimited"
	echo " sudo swapoff -a; sudo swapon -a : clean swap"
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
