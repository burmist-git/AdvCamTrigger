#!/bin/bash

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d            : default"
    echo " [0] --conv_to_csv : convert to csv"
    echo " [1]   part_type   : (g/gd/e/p)"
    echo " [0] --size        : get data size"
    echo " [0] --ls          : ls data"
    echo " [0] -h            : print help"
}

function conv_to_csv {
    #
    #corsika_run1.header.pkl_0
    #corsika_run1.pe_info.pkl_0
    #../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/pkl/0000/
    #
    particle=$1
    particlein=$2
    folderID=$3
    pklpath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/pkl/$folderID/"
    csvpath="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/$folderID/"
    mkdir -p $csvpath
    #
    echo "particle   $particle"
    echo "particlein $particlein"
    echo "folderID   $folderID"
    echo "pklpath    $pklpath"
    echo "csvpath    $csvpath"
    #
    for filename_header_pkl in `ls -lrt $pklpath/corsika_run*.header.pkl_* | awk {'print $9'}`
    do
	filename_header_csv=$csvpath`basename $filename_header_pkl`'.csv'
	python pkl_tocsv.py header $filename_header_pkl $filename_header_csv
    done
    #
    for filename_pe_pkl in `ls -lrt $pklpath/corsika_run*.pe_info.pkl_* | awk {'print $9'}`
    do
	filename_pe_csv=$csvpath`basename $filename_pe_pkl`'.csv'
	python pkl_tocsv.py pe $filename_pe_pkl $filename_pe_csv
    done
}

if [ $# -eq 0 ]; then
    printHelp
else
    if [ "$1" = "--conv_to_csv" ]; then
	if [ $# -eq 2 ]; then
	    if [ "$2" = "g" ]; then
		conda activate pyeventio
		particle="gamma"
		particlein="gamma_on_nsb_1x"
		i_start=0
		i_stop=0
		#
		for i in `seq $i_start $i_stop`
		do
		    folderID=`printf %04d $i`
		    #echo "$folderID"
		    conv_to_csv $particle $particlein $folderID
		done
	    elif [ "$2" = "gd" ]; then
		conda activate pyeventio
		particle="gamma_diffuse"
		particlein="gamma_diffuse_nsb_1x"
		i_start=0
		i_stop=0
		#
		for i in `seq $i_start $i_stop`
		do
		    folderID=`printf %04d $i`
		    conv_to_csv $particle $particlein $folderID
		done
	    elif [ "$2" = "e" ]; then
		conda activate pyeventio
		particle="electron"
		particlein="electron_nsb_1x"
		i_start=0
		i_stop=0
		#
		for i in `seq $i_start $i_stop`
		do
		    folderID=`printf %04d $i`
		    conv_to_csv $particle $particlein $folderID
		done
	    elif [ "$2" = "p" ]; then
		conda activate pyeventio
		particle="proton"
		particlein="proton_nsb_1x"
		i_start=0
		i_stop=0
		#
		for i in `seq $i_start $i_stop`
		do
		    folderID=`printf %04d $i`
		    conv_to_csv $particle $particlein $folderID
		done
	    else
		printHelp
	    fi
	else
            printHelp
	fi
    elif [ "$1" = "--size" ]; then
	#
	rootPathsimteldata="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/"
	#
	particlein="gamma_on_nsb_1x"
	csvpath="$rootPathsimteldata$particlein/csv/"
	du -hs $csvpath
	#
	particlein="gamma_diffuse_nsb_1x"
	csvpath="$rootPathsimteldata$particlein/csv/"
	du -hs $csvpath
	#
	particlein="electron_nsb_1x"
	csvpath="$rootPathsimteldata$particlein/csv/"
	du -hs $csvpath
	#
	particlein="proton_nsb_1x"
	csvpath="$rootPathsimteldata$particlein/csv/"
	du -hs $csvpath
    elif [ "$1" = "--ls" ]; then
	#
	rootPathsimteldata="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/"
	#
	particlein="gamma_on_nsb_1x"
	csvpath="$rootPathsimteldata$particlein/csv/"
	ls $csvpath
	#
	particlein="gamma_diffuse_nsb_1x"
	csvpath="$rootPathsimteldata$particlein/csv/"
	ls $csvpath
	#
	particlein="electron_nsb_1x"
	csvpath="$rootPathsimteldata$particlein/csv/"
	ls $csvpath
	#
	particlein="proton_nsb_1x"
	csvpath="$rootPathsimteldata$particlein/csv/"
	ls $csvpath
    elif [ "$1" = "-d" ]; then
        printHelp
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
