#!/bin/bash
#SBATCH --job-name crgen%j
#SBATCH --error /home/users/b/burmistr//job_error/crgen_%j.error
#SBATCH --output /home/users/b/burmistr//job_output/output_%j.output
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition public-cpu
#SBATCH --time 0-03:00:00

#source /home/users/b/burmistr/terzina_photon_propagation/setupEnv.sh -d

################################################
#  E                                           #
################################################
#  0                    1              1.58489 #
#  1              1.58489              2.51189 #
#  2              2.51189              3.98107 #
#  3              3.98107              6.30957 #
#  4              6.30957                   10 #
#  5                   10              15.8489 #
#  6              15.8489              25.1189 #
#  7              25.1189              39.8107 #
#  8              39.8107              63.0957 #
#  9              63.0957                  100 #
# 10                  100              158.489 #
# 11              158.489              251.189 #
# 12              251.189              398.107 #
# 13              398.107              630.957 #
# 14              630.957                 1000 #
# 15                 1000              1584.89 #
# 16              1584.89              2511.89 #
# 17              2511.89              3981.07 #
# 18              3981.07              6309.57 #
# 19              6309.57                10000 #
# 20                10000              15848.9 #
# 21              15848.9              25118.9 #
# 22              25118.9              39810.7 #
# 23              39810.7              63095.7 #
# 24              63095.7               100000 #
################################################
# theta                                        #
################################################
#  0                    0                    1 #
#  1                    1                    2 #
#  2                    2                    3 #
#  3                    3                    4 #
#  4                    4                    5 #
#  5                    5                    6 #
#  6                    6                    7 #
#  7                    7                    8 #
#  8                    8                    9 #
#  9                    9                   10 #
################################################
#  r                                           #
################################################
#  0                    0                  150 #
#  1                  150                  250 #
#  2                  250                  350 #
#  3                  350                  500 #
#  4                  500                  650 #
#  5                  650                  800 #
#  6                  800                 1000 #
#  7                 1000                 1300 #
#  8                 1300                 1600 #
#  9                 1600                 2000 #
################################################



function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d     : default"
    echo " [1]        : particle type (g,gd,e,p)"
    echo " [2]        : Ebin   - [0:24]"
    echo " [3]        : Thbin  - [0:9]"
    echo " [4]        : rbin   - [0:9]"
    echo " [5]        : jobID  - [0:199]"
    echo " [0] -NGB   : NGB"
    echo " [1]        : jobID  - [0:199] (ex: 0000  0001 ...)"
    echo " [0] -c     : recompile"
    echo " [0] -h     : print help"
}

if [ $# -eq 0 ]; then
    printHelp
else
    if [ "$1" = "-d" ]; then
	if [ $# -eq 6 ]; then
	    if [ "$2" = "g" ]; then
		inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/root/"
		outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/trgA/"
	    elif [ "$2" = "gd" ]; then
		inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_diffuse_nsb_1x/root/"
		outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_diffuse_nsb_1x/trgA/"
	    elif [ "$2" = "e" ]; then
		inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/electron_nsb_1x/root/"
		outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/electron_nsb_1x/trgA/"
	    elif [ "$2" = "p" ]; then
		inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/"
		outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/trgA/"
	    fi
	    binE=$3
	    binTheta=$4
	    binDist=$5
	    jobID=$6
	    inRootFile=$inRootFilePref$jobID"/corsika_"$jobID"ID.root"
	    mkdir -p $outHistFPref$jobID
	    outHistF=$outHistFPref$jobID"/hist_trgA_corsika_"$binE"binE_"$binTheta"binTheta_"$binDist"binDist_"$jobID"ID.root"
	    npe_min=20
	    npe_max=200
	    nEv_max=100
	    rndseed=`date +%N`
	    echo "inRootFile $inRootFile"	
	    echo "outHistF   $outHistF"
	    echo "binE       $binE"
	    echo "binTheta   $binTheta"
	    echo "binDist    $binDist"
	    echo "jobID      $jobID"
	    echo "npe_min    $npe_min"	
	    echo "npe_max    $npe_max"
	    echo "nEv_max    $nEv_max"
	    echo "rndseed    $rndseed"
	    #live
	    #./runana 111 $inRootFile $outHistF $binE $binTheta $binDist $npe_min $npe_max $nEv_max $rndseed
	    #screen
	    screenName='tr'$jobID
            echo "$screenName"
            screen -S $screenName -L -d -m ./runana 111 $inRootFile $outHistF $binE $binTheta $binDist $npe_min $npe_max $nEv_max $rndseed
	    #srun
	else
	    printHelp	    
	fi	
    elif [ "$1" = "-NGB" ]; then
	if [ $# -eq 2 ]; then
	    inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/"
	    outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgA/"
	    jobID=$2
	    inRootFile=$inRootFilePref$jobID"/corsika_"$jobID"ID.root"
	    outHistF=$outHistFPref"/hist_trgA_corsika_"$jobID"ID.root"
	    nEv_max=1000
	    rndseed=`date +%N`
	    echo "inRootFile $inRootFile"	
	    echo "outHistF   $outHistF"
	    echo "jobID      $jobID"
	    echo "nEv_max    $nEv_max"
	    echo "rndseed    $rndseed"
	    #live
	    #./runana 112 $inRootFile $outHistF $nEv_max $rndseed
	    #screen
	    screenName='tr'$jobID
            echo "$screenName"
            screen -S $screenName -L -d -m ./runana 112 $inRootFile $outHistF $nEv_max $rndseed
	    #srun
	else
	    printHelp
	fi
    elif [ "$1" = "-c" ]; then
	make clean; make -f Makefileana clean ; make -j; make -f Makefileana -j		
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi


