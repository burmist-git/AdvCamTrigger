#!/bin/bash

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -m      : muons"
    echo " [0] --mfast : muons fast"
    echo " [0] -c      : recompile"
    echo " [0] --vis   : vis"
    echo " [0] -h      : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-m" ]; then
	inRootFile="../scratch/simtel_data/muon/root/run1_muon.root"
	outHistF="../scratch/simtel_data/muon/hist/hist_run1_muon.root"
	#outHistF="../scratch/simtel_data/muon/hist/hist_run1_muon_evID_cut.root"
	./runanamuon 1 $inRootFile $outHistF
    elif [ "$1" = "--mfast" ]; then
	inRootFile="../scratch/simtel_data/muon/root/run1_muon.root"
	outHistF="../scratch/simtel_data/muon/hist/hist_fast_run1_muon.root"
	./runanamuon 11 $inRootFile $outHistF
    elif [ "$1" = "--vis" ]; then
	inRootFile="../scratch/simtel_data/muon/root/run1_muon.root"
	outHistF="../scratch/simtel_data/muon/hist/hist_vis_run1_muon.root"
      	#evID=450
	#evID=403
	#evID=885
	#evID=726
	#evID=23
	#evID=3018
	#evID=1950
	#evID=2066
	#evID=2756
	#evID=24
	#evID=465
	#evID=225
	evID=5371
	particle_type_name="muon"
	./runanamuon 2 $inRootFile $outHistF $evID $particle_type_name
	evID=2039
	particle_type_name="muon"
	./runanamuon 2 $inRootFile $outHistF $evID $particle_type_name
    elif [ "$1" = "-c" ]; then
	make clean; make -f Makefileana clean ; make -j; make -f Makefileana -j		
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak    "I have done"


