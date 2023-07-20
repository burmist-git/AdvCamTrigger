#!/bin/bash

#'nightsky_background=all:0.001'
#inRootFile="../compressed_data/no_nsb_cut/gamma/corsika_run307.compressed.root"
#outHistF="./hist_corsika_run307_gamma_no_nsb_cut.root"
#inRootFile="../compressed_data/no_nsb_cut/proton/corsika_run307.compressed.root"
#outHistF="./hist_corsika_run307_proton_no_nsb_cut.root"
#'nightsky_background=all:0.386'
inRootFile="../compressed_data/gamma/corsika_run307.compressed.root"
outHistF="./hist_corsika_run307.root"

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
    echo " [0] -th : test cam sipmCameraHist"
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
	#./runana 4 $inRootFile $outHistF 896 gamma    # 162 pe
	#./runana 4 $inRootFile $outHistF 929 gamma    # 220 pe
	#./runana 4 $inRootFile $outHistF 3140 gamma   # 180 pe
	#./runana 4 $inRootFile $outHistF 3139 gamma   # 210 pe
	./runana 4 $inRootFile $outHistF 1642 gamma   # 102 pe
	./runana 4 $inRootFile $outHistF 1673 gamma   # 98 pe
	./runana 4 $inRootFile $outHistF 1699 gamma   # 100 pe
	./runana 4 $inRootFile $outHistF 1708 gamma   # 101 pe
	#./runana 4 $inRootFile $outHistF 375 gamma   # 43 pe
	#./runana 4 $inRootFile $outHistF 397 gamma   # 43 pe
	#./runana 4 $inRootFile $outHistF 3058 gamma  # 43 pe
	#./runana 4 $inRootFile $outHistF 2153 gamma  # 43 pe
	#./runana 4 $inRootFile $outHistF 3929 gamma  # 42 pe
	#./runana 4 $inRootFile $outHistF 3932 gamma  # 44 pe
	#./runana 4 $inRootFile $outHistF 4067 gamma   # 42 pe
	#./runana 4 $inRootFile $outHistF 4064 gamma   # 40 pe
	#./runana 4 $inRootFile $outHistF 1222 gamma   # 27 pe
	#./runana 4 $inRootFile $outHistF 1351 gamma    # 32 pe
	#./runana 4 $inRootFile $outHistF 1374 gamma    # 30 pe
	#./runana 4 $inRootFile $outHistF 2993 gamma    # 29 pe
	#./runana 4 $inRootFile $outHistF 3792 gamma    # 30 pe
	#./runana 4 $inRootFile $outHistF 3799 gamma    # 30 pe
	#./runana 4 $inRootFile $outHistF 4540 gamma    # 30 pe
	#./runana 4 $inRootFile $outHistF 4 gamma    # 20 pe
	#./runana 4 $inRootFile $outHistF 111 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 118 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 142 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 405 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 872 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 848 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 952 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 959 gamma # 20 pe
    elif [ "$1" = "-th" ]; then
	./runana 5
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak "I have done"
