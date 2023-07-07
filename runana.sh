#!/bin/bash

#'nightsky_background=all:0.001'
inRootFile="../compressed_data/no_nsb_cut/gamma/corsika_run307.compressed.root"
outHistF="./hist_corsika_run307_gamma_no_nsb_cut.root"
#inRootFile="../compressed_data/no_nsb_cut/proton/corsika_run307.compressed.root"
#outHistF="./hist_corsika_run307_proton_no_nsb_cut.root"
#'nightsky_background=all:0.386'
#inRootFile="../compressed_data/gamma/corsika_run307.compressed.root"
#outHistF="./hist_corsika_run307.root"

#short
#'nightsky_background=all:0.001'
#inRootFile="../compressed_data/no_nsb_cut/gamma/corsika_run307.compressed.short.root"
#outHistF="./hist_short_corsika_run307_no_nsb_cut.root"
#'nightsky_background=all:0.386'
#inRootFile="../compressed_data/gamma/corsika_run307.compressed.root"
#outHistF="./hist_corsika_run307.root"

make -f Makefileana clean; make -f Makefileana runana;

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d  : single root file"
    echo " [0] -ds : single root file (short format)"
    echo " [0] -tw : test waveforms"
    echo " [0] -h  : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	./runana 1 $inRootFile $outHistF
    elif [ "$1" = "-ds" ]; then
	./runana 3 $inRootFile $outHistF
    elif [ "$1" = "-tw" ]; then
	./runana 4 $inRootFile $outHistF 2500
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak "I have done"
