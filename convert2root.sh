#!/bin/bash

# Thu 11 May 16:09:40 CEST 2023
# Autor: Leonid Burmistrov

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] --convert_gm_test  : convert gammas         (test)"
    echo " [0] --convert_gd_test  : convert gammas diffuse (test)"
    echo " [0] --convert_el_test  : convert electrons      (test)"
    echo " [0] --convert_pr_test  : convert protons        (test)"
    echo " [0] --load_modules     : load_modules"
    echo " [0] --convert_gm       : convert gammas"
    echo " [0] --convert_gd       : convert gammas diffuse"
    echo " [0] --convert_el       : convert electrons"
    echo " [0] --convert_pr       : convert protons"
    echo " [0] --inc_reduce_stack : increase reduce stack memory"
    echo " [0] -h                 : print help"
}

function creat_csv2toot_list {
    #corsika_run1.header.pkl_0.csv
    #corsika_run1.pe_info.pkl_0.csv
    csv_dir=$1
    n_run_header=`ls $csv_dir/corsika_run*.header.pkl_0.csv | wc -w`
    n_run_pe=`ls $csv_dir/corsika_run*.pe_info.pkl_0.csv | wc -w`
    list="$csv_dir/list"
    rm -rf $list
    echo " creat_csv2toot_list     --> $csv_dir"
    echo " n_run_header            --> $n_run_header"
    echo " n_run_pe                --> $n_run_pe"
    echo " list                    --> $list"
    #
    for run_i in `seq 1 $n_run_header`
    do
	n_files_header=`ls $csv_dir/corsika_run$run_i.header.pkl_*.csv | wc -w`
	n_files_pe=`ls $csv_dir/corsika_run$run_i.pe_info.pkl_*.csv | wc -w`
	echo " run: $run_i n_files_header   --> $n_files_header"
	echo "         n_run_pe         --> $n_files_pe"
	for pkl_i in `seq 0 $(echo "$n_run_header-1" | bc )`
	do
	    header_csv_file="corsika_run$run_i.header.pkl_$pkl_i.csv"
	    pe_info_csv_file="corsika_run$run_i.pe_info.pkl_$pkl_i.csv"	    
	    if [ -f "$csv_dir/$header_csv_file" ]; then
		if [ -f "$csv_dir/$pe_info_csv_file" ]; then
		    #echo "$header_csv_file $pe_info_csv_file" | tee -a $csv_dir/list
		    echo "$header_csv_file $pe_info_csv_file" >> $list
		fi
		#echo "$header_csv_file" | tee -a $list
		#echo "$header_csv_file" >> $list
	    fi
	done
    done
}

function load_modules {
    module load GCC/7.3.0-2.30 GCCcore/7.3.0 OpenMPI/3.1.1 ROOT/6.14.06-Python-2.7.15 Geant4/10.5 CMake/3.11.4
}
    
if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "--convert_gm_test" ]; then
	make clean; make;
	particle="gamma"
	particlein="gamma_on_nsb_1x"
	#
	#header_file="../compressed_data/no_nsb_cut/$particle/corsika_run307.header.csv"
	#pe_info_file="../compressed_data/no_nsb_cut/$particle/corsika_run307.pe_info.csv"
	#wf_file_list="./wf_file_list_no_nsb_cut.list"
	#outputRootFile="../compressed_data/no_nsb_cut/$particle/corsika_run307.compressed.root"
	#outputRootFile="../compressed_data/no_nsb_cut/$particle/corsika_run307.compressed.short.root"
	#
	#header_file="../compressed_data/$particle/corsika_run307.header.csv"
	#pe_info_file="../compressed_data/$particle/corsika_run307.pe_info.csv"
	#wf_file_list="./wf_file_list.list"
	#outputRootFile="../compressed_data/$particle/corsika_run307.compressed.root"
	#
	#./convert2root 0 $header_file $pe_info_file $wf_file_list $outputRootFile
	#
	header_file="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/0000/corsika_run1.header.pkl_0.csv"
	pe_info_file="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/0000/corsika_run1.pe_info.pkl_0.csv"
	outputRootFilePath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/root/0000/"
	mkdir -p $outputRootFilePath
	outputRootFile="$outputRootFilePath/corsika_0000ID.root"
	./convert2root 1 $header_file $pe_info_file $outputRootFile
    elif [ "$1" = "--convert_gd_test" ]; then
	make clean; make;
	particle="gamma_diffuse"
	particlein="gamma_diffuse_nsb_1x"
	header_file="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/0000/corsika_run1.header.pkl_0.csv"
	pe_info_file="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/0000/corsika_run1.pe_info.pkl_0.csv"
	outputRootFilePath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/root/0000/"
	mkdir -p $outputRootFilePath
	outputRootFile="$outputRootFilePath/corsika_0000ID.root"
	./convert2root 1 $header_file $pe_info_file $outputRootFile
    elif [ "$1" = "--convert_el_test" ]; then
	make clean; make;
	particle="electron"
	particlein="electron_nsb_1x"
	header_file="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/0000/corsika_run1.header.pkl_0.csv"
	pe_info_file="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/0000/corsika_run1.pe_info.pkl_0.csv"
	outputRootFilePath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/root/0000/"
	mkdir -p $outputRootFilePath
	outputRootFile="$outputRootFilePath/corsika_0000ID.root"
	./convert2root 1 $header_file $pe_info_file $outputRootFile
    elif [ "$1" = "--convert_pr_test" ]; then
	make clean; make;
	particle="proton"
	particlein="proton_nsb_1x"
	header_file="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/0000/corsika_run1.header.pkl_0.csv"
	pe_info_file="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/0000/corsika_run1.pe_info.pkl_0.csv"
	outputRootFilePath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/root/0000/"
	mkdir -p $outputRootFilePath
	outputRootFile="$outputRootFilePath/corsika_0000ID.root"
	./convert2root 1 $header_file $pe_info_file $outputRootFile
    elif [ "$1" = "--load_modules" ]; then	
	load_modules
    elif [ "$1" = "--convert_gm" ]; then
	load_modules
	#
	particle="gamma"
	particlein="gamma_on_nsb_1x"
	i_start=0
	i_stop=0
	#
	for i in `seq $i_start $i_stop`
	do
	    folderID=`printf %04d $i`
	    csv_dir="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/$folderID/"
	    outputRootFilePath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/root/$folderID/"
	    mkdir -p $outputRootFilePath
	    creat_csv2toot_list $csv_dir
	    ./convert2root 2 $csv_dir "$csv_dir/list" $outputRootFilePath/"corsika_"$folderID"ID.root"
	done
    elif [ "$1" = "--convert_gd" ]; then
	load_modules
	#
	particle="gamma_diffuse"
	particlein="gamma_diffuse_nsb_1x"
	i_start=0
	i_stop=0
	#
	for i in `seq $i_start $i_stop`
	do
	    folderID=`printf %04d $i`
	    csv_dir="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/$folderID/"
	    outputRootFilePath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/root/$folderID/"
	    mkdir -p $outputRootFilePath
	    creat_csv2toot_list $csv_dir
	    ./convert2root 2 $csv_dir "$csv_dir/list" $outputRootFilePath/"corsika_"$folderID"ID.root"
	done
    elif [ "$1" = "--convert_el" ]; then
	load_modules
	#
	particle="electron"
	particlein="electron_nsb_1x"
	i_start=0
	i_stop=0
	#
	for i in `seq $i_start $i_stop`
	do
	    folderID=`printf %04d $i`
	    csv_dir="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/$folderID/"
	    outputRootFilePath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/root/$folderID/"
	    mkdir -p $outputRootFilePath
	    creat_csv2toot_list $csv_dir
	    ./convert2root 2 $csv_dir "$csv_dir/list" $outputRootFilePath/"corsika_"$folderID"ID.root"
	done
    elif [ "$1" = "--convert_pr" ]; then
	load_modules
	#
	particle="proton"
	particlein="proton_nsb_1x"
	i_start=0
	i_stop=0
	#
	for i in `seq $i_start $i_stop`
	do
	    folderID=`printf %04d $i`
	    csv_dir="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/$folderID/"
	    outputRootFilePath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/root/$folderID/"
	    mkdir -p $outputRootFilePath
	    creat_csv2toot_list $csv_dir
	    ./convert2root 2 $csv_dir "$csv_dir/list" $outputRootFilePath/"corsika_"$folderID"ID.root"
	done
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
