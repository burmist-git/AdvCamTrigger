#!/bin/bash

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d : root files list"
    echo " [0] -t : single root file test"
    echo " [0] -c : recompile"
    echo " [0] -h : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	fileWithListOfRootFiles="proton_st.list"
	outRootFileHist="hist_proton_st.root"
	./runanastereo 0 $fileWithListOfRootFiles $outRootFileHist		
    elif [ "$1" = "-t" ]; then
	inRootFile="corsika_proton_run1.root"
	outRootFileHist="hist_corsika_proton_run1.root"
	./runanastereo 1 $inRootFile $outRootFileHist	
    elif [ "$1" = "-c" ]; then
	make -f Makefileana cleanstereo; make -f Makefileana runanastereo  -j
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak    "I have done"


