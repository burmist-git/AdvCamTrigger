#!/bin/bash

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d            : default"
    echo " [0] --conv_to_pkl : convert to pkl"
    echo " [1]   part_type   : (g/gd/e/p)"
    echo " [2]   i_start     : i_start"
    echo " [3]   i_stop      : i_stop"
    echo " [0] --size        : get data size"
    echo " [0] --ls          : ls data"
    echo " [0] --conv_all_gamma         : convert to pkl all gamma"
    echo " [0] --conv_all_gamma_diffuse : convert to pkl all gamma_diffuse"
    echo " [0] --conv_all_electron      : convert to pkl all electron"
    echo " [0] --conv_all_proton        : convert to pkl all proton"
    echo " [0] --rm_900_gamma_simtel_gz          : rm 900 gamma simtel.gz files"
    echo " [0] --rm_1800_gamma_diffuse_simtel_gz : rm 1800 gamma diffuse simtel.gz files"
    echo " [0] --rm_1800_electron_simtel_gz      : rm 1800 electron simtel.gz files"
    echo " [0] --rm_4500_proton_simtel_gz        : rm 4500 proton simtel.gz files"
    echo " [0] -h                                : print help"
}

function rm_simtel_gz {
    particle=$1
    particlein=$2
    n_rm_files=$3
    n_files_max=$4
    simtelpath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/output/"
    echo "particle   $particle"
    echo "particlein $particlein"
    echo "n_rm_files $n_rm_files"

    echo "simtelpath $simtelpath"
    python gen_rand_tuple.py $n_rm_files $n_files_max > "ID_to_rm."$particle
    #while read -r runID;
    #do rm $simtelpath/corsika_run$runID.simtel.gz;
    #done < "ID_to_rm."$particle
}

function rm_percentage_simtel_gz {
    particle=$1
    particlein=$2
    percentage=$3
    simtelpath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/output/"
    echo "particle   $particle"
    echo "particlein $particlein"
    echo "percentage $percentage"
    echo "simtelpath $simtelpath"
    #
    for f in $simtelpath/*
    do
	if [ $((1 + $RANDOM % 100)) \> $percentage ];
	then
	    echo "$f"
	else
	    rm $f
	fi;
    done
}

function conv_to_pkl {
    #
    #corsika_run1.simtel.gz
    #corsika_run1.header.pkl
    #corsika_run1.pe_info.pkl
    #
    particle=$1
    particlein=$2
    ind_start=$3
    ind_stop=$4
    #simtelpath="../scratch/data_nagaia/data/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/output/"
    #simtelpath="../simtel_data/$particle/data/"
    simtelpath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/output/"
    pklpath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/pkl/"
    mkdir -p $pklpath
    #
    echo "particle   $particle"
    echo "particlein $particlein"
    echo "ind_start  $ind_start"
    echo "ind_stop   $ind_stop"
    echo "simtelpath $simtelpath"
    #
    #ind_stop=$4
    for i in `seq $ind_start $ind_stop`
    do
	filename=$simtelpath"corsika_run$i.simtel.gz"
	datafilein=$filename
	#
	folderID=`echo "($i - 1) / 100" | bc`
	folderName=$pklpath`printf %04d $folderID`
	mkdir -p $folderName
	echo "folderName = $folderName"
	#
        headerout="$folderName/corsika_run$i.header.pkl"
        pe_info_out="$folderName/corsika_run$i.pe_info.pkl"
	python pkl_pyeventio.py $datafilein $headerout $pe_info_out
    done
}

if [ $# -eq 0 ]; then
    printHelp
else
    if [ "$1" = "--conv_to_pkl" ]; then
	if [ $# -eq 4 ]; then
	    if [ "$2" = "g" ]; then
		particle="gamma"
		particlein="gamma_on_nsb_1x"
		i_start=$3
		i_stop=$4
		conv_to_pkl $particle $particlein $i_start $i_stop
	    elif [ "$2" = "gd" ]; then
		particle="gamma_diffuse"
		particlein="gamma_diffuse_nsb_1x"
		i_start=$3
		i_stop=$4
		conv_to_pkl $particle $particlein $i_start $i_stop
	    elif [ "$2" = "e" ]; then
		particle="electron"
		particlein="electron_nsb_1x"
		i_start=$3
		i_stop=$4
		conv_to_pkl $particle $particlein $i_start $i_stop
	    elif [ "$2" = "p" ]; then
		particle="proton"
		particlein="proton_nsb_1x"
		i_start=$3
		i_stop=$4
		conv_to_pkl $particle $particlein $i_start $i_stop
	    else
		printHelp
	    fi
	else
            printHelp
	fi
    elif [ "$1" = "--conv_all_gamma" ]; then
	conda activate pyeventio
	particle="gamma"
	particlein="gamma_on_nsb_1x"
	i_start=1
	i_stop=1000
	conv_to_pkl $particle $particlein $i_start $i_stop
    elif [ "$1" = "--conv_all_gamma_diffuse" ]; then
	conda activate pyeventio
	particle="gamma_diffuse"
	particlein="gamma_diffuse_nsb_1x"
	i_start=1901
	i_stop=2000
	conv_to_pkl $particle $particlein $i_start $i_stop
    elif [ "$1" = "--conv_all_electron" ]; then
	conda activate pyeventio
	particle="electron"
	particlein="electron_nsb_1x"
	i_start=1
	i_stop=2000
	conv_to_pkl $particle $particlein $i_start $i_stop
    elif [ "$1" = "--conv_all_proton" ]; then
	conda activate pyeventio
	particle="proton"
	particlein="proton_nsb_1x"
	i_start=1
	i_stop=5000
	conv_to_pkl $particle $particlein $i_start $i_stop	
    elif [ "$1" = "--rm_900_gamma_simtel_gz" ]; then
	particle="gamma"
	particlein="gamma_on_nsb_1x"
	n_rm_files=900
	n_files_max=1000
	percentage=47
	#rm_simtel_gz $particle $particlein $n_rm_files $n_files_max
	#rm_percentage_simtel_gz $particle $particlein $percentage
    elif [ "$1" = "--rm_1800_gamma_diffuse_simtel_gz" ]; then
	particle="gamma_diffuse"
	particlein="gamma_diffuse_nsb_1x"
	n_rm_files=1800
	n_files_max=2000
	percentage=83
	#rm_simtel_gz $particle $particlein $n_rm_files $n_files_max
	#rm_percentage_simtel_gz $particle $particlein $percentage
    elif [ "$1" = "--rm_1800_electron_simtel_gz" ]; then
	particle="electron"
	particlein="electron_nsb_1x"
	n_rm_files=1800
	n_files_max=2000
	percentage=61
	#rm_simtel_gz $particle $particlein $n_rm_files $n_files_max
	#rm_percentage_simtel_gz $particle $particlein $percentage
    elif [ "$1" = "--rm_4500_proton_simtel_gz" ]; then
	particle="proton"
	particlein="proton_nsb_1x"
	n_rm_files=4500
	n_files_max=5000
	percentage=80
	#rm_simtel_gz $particle $particlein $n_rm_files $n_files_max
	#rm_percentage_simtel_gz $particle $particlein $percentage
    elif [ "$1" = "--size" ]; then
	#
	rootPathsimteldata="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/"
	#
	particlein="gamma_on_nsb_1x"
	pklpath="$rootPathsimteldata$particlein/pkl/"
	du -hs $pklpath
	#
	particlein="gamma_diffuse_nsb_1x"
	pklpath="$rootPathsimteldata$particlein/pkl/"
	du -hs $pklpath
	#
	particlein="electron_nsb_1x"
	pklpath="$rootPathsimteldata$particlein/pkl/"
	du -hs $pklpath
	#
	particlein="proton_nsb_1x"
	pklpath="$rootPathsimteldata$particlein/pkl/"
	du -hs $pklpath
    elif [ "$1" = "--ls" ]; then
	#
	rootPathsimteldata="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/"
	#
	particlein="gamma_on_nsb_1x"
	pklpath="$rootPathsimteldata$particlein/pkl/"
	ls $pklpath
	#
	particlein="gamma_diffuse_nsb_1x"
	pklpath="$rootPathsimteldata$particlein/pkl/"
	ls $pklpath
	#
	particlein="electron_nsb_1x"
	pklpath="$rootPathsimteldata$particlein/pkl/"
	ls $pklpath
	#
	particlein="proton_nsb_1x"
	pklpath="$rootPathsimteldata$particlein/pkl/"
	ls $pklpath
    elif [ "$1" = "-d" ]; then
        printHelp
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
