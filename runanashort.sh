#!/bin/bash

#'nightsky_background=all:0.001'
inRootFile="../compressed_data/no_nsb_cut/gamma/corsika_run307.compressed.short.root"
outHistF="./hist_short_corsika_run307_no_nsb_cut.root"

#'nightsky_background=all:0.386'
#inRootFile="../compressed_data/gamma/corsika_run307.compressed.root"
#outHistF="./hist_corsika_run307.root"

make -f Makefileana clean; make -f Makefileana runanashort;

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d  : single root file"
    echo " [0] -h  : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	./runanashort 1 $inRootFile $outHistF
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak "I have done"
