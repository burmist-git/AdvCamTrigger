#!/bin/bash
#SBATCH --job-name crgen%j
#SBATCH --error /srv/beegfs/scratch/users/b/burmistr/pyeventio_example/job_error/crgen_%j.error
#SBATCH --output /srv/beegfs/scratch/users/b/burmistr/pyeventio_example/job_output/output_%j.output
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition public-cpu
#SBATCH --time 0-03:00:00

#simHomeDir="/home/users/b/burmistr/pyeventio_example/"
simHomeDir="./"

source $simHomeDir/setupEnv.sh -d
source $simHomeDir/convert2root.sh --set_unlimited_mem

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
    echo " [0] -d                    : default"
    echo " [1]                       : particle type (g,gd,e,p)"
    echo " [2]                       : Ebin          - [0:24]"
    echo " [3]                       : Thbin         - [0:9]"
    echo " [4]                       : rbin          - [0:9]"
    echo " [5]                       : jobID         - [0:199]"
    echo " [6]                       : data_chunk_ID - [0:19]"
    echo " [7]                       : jobType       - [sbatch:screen:live]"
    echo " [0] -NGB                  : NGB"
    echo " [1]                       : jobID - [0:199] (ex: 0000  0001 ...)"
    echo " [2]                       : jobType       - [sbatch:screen:live]"
    echo " [0] -test_live_NSB        : test live (NSB)"
    echo " [1]                       : nEv"
    echo " [0] -test_live_NSB_k_dist : test live (NSB) k-dist plot"
    echo " [0] -test_live            : test live"
    echo " [1]                       : nEv"
    echo " [0] -test_srun            : test srun"
    echo " [1]                       : nEv"
    echo " [0] -c                    : recompile"
    echo " [0] -h                    : print help"
}

if [ $# -eq 0 ]; then
    printHelp
else
    if [ "$1" = "-d" ]; then
	if [ $# -eq 8 ]; then
	    if [ "$2" = "g" ]; then
		inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/root/"
		#outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/trgA/"
		outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/trgA_test/"
		rsimulation=800
	    elif [ "$2" = "gd" ]; then
		inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_diffuse_nsb_1x/root/"
		#outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_diffuse_nsb_1x/trgA/"
		outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_diffuse_nsb_1x/trgA_test/"
		rsimulation=1000
	    elif [ "$2" = "e" ]; then
		inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/electron_nsb_1x/root/"
		#outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/electron_nsb_1x/trgA/"
		outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/electron_nsb_1x/trgA_test/"
		rsimulation=1000
	    elif [ "$2" = "p" ]; then
		inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/"
		#outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/trgA/"
		outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/trgA_test/"
		rsimulation=1500
	    fi
	    binE=$3
	    binTheta=$4
	    binDist=$5
	    jobID=$6
	    data_chunk_ID=$7
	    typeJob=$8
	    #
	    inRootFile=$inRootFilePref$jobID"/corsika_"$jobID"ID.root"
	    mkdir -p $outHistFPref$jobID
	    outHistF=$outHistFPref$jobID"/hist_trgA_corsika_"$binE"binE_"$binTheta"binTheta_"$binDist"binDist_"$jobID"ID.root"
	    npe_min=20
	    npe_max=10000
	    #nEv_max=15000
	    nEv_max=10
	    rndseed=`date +%N`
	    trgSetup="trg_setup.conf"
	    echo "inRootFile    $inRootFile"
	    echo "outHistF      $outHistF"
	    echo "binE          $binE"
	    echo "binTheta      $binTheta"
	    echo "binDist       $binDist"
	    echo "jobID         $jobID"
	    echo "npe_min       $npe_min"	
	    echo "npe_max       $npe_max"
	    echo "nEv_max       $nEv_max"
	    echo "rndseed       $rndseed"
	    echo "data_chunk_ID $data_chunk_ID"
	    echo "rsimulation   $rsimulation"
	    echo "trgSetup      $trgSetup"	    
	    #
	    #live
	    #./runana 111 $inRootFile $outHistF $binE $binTheta $binDist $npe_min $npe_max $nEv_max $rndseed $data_chunk_ID
	    #screen
	    #screenName='tr'$jobID
            #echo "$screenName"
            #screen -S $screenName -L -d -m ./runana 111 $inRootFile $outHistF $binE $binTheta $binDist $npe_min $npe_max $nEv_max $rndseed $data_chunk_ID
	    #srun
	    #
	    if [ "$typeJob" = "sbatch" ]; then
		srun $simHomeDir/runana 111 $inRootFile $outHistF $binE $binTheta $binDist $npe_min $npe_max $nEv_max $rndseed $data_chunk_ID $rsimulation $trgSetup
	    elif [ "$typeJob" = "screen" ]; then
		screenName='tr'$jobID
		echo "$screenName"
		screen -S $screenName -L -d -m $simHomeDir/runana 111 $inRootFile $outHistF $binE $binTheta $binDist $npe_min $npe_max $nEv_max $rndseed $data_chunk_ID $rsimulation $trgSetup
	    elif [ "$typeJob" = "live" ]; then
		$simHomeDir/runana 111 $inRootFile $outHistF $binE $binTheta $binDist $npe_min $npe_max $nEv_max $rndseed $data_chunk_ID $rsimulation $trgSetup
	    else
		printHelp
	    fi
	else
	    printHelp	    
	fi	
    elif [ "$1" = "-NGB" ]; then
	if [ $# -eq 3 ]; then
	    inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/"
	    #outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgA/"
	    #outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_386MHz/trgA/"
	    outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgA_test/"
	    #outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_386MHz/trgA_test/"
	    mkdir -p $outHistFPref
	    jobID=$2
	    typeJob=$3
	    inRootFile=$inRootFilePref$jobID"/corsika_"$jobID"ID.root"
	    outHistF=$outHistFPref"/hist_trgA_corsika_"$jobID"ID.root"
	    nEv_max=10
	    rndseed=`date +%N`
	    trgSetup="trg_setup.conf"
	    echo "inRootFile $inRootFile"	
	    echo "outHistF   $outHistF"
	    echo "jobID      $jobID"
	    echo "nEv_max    $nEv_max"
	    echo "rndseed    $rndseed"
	    echo "typeJob    $typeJob"
	    echo "trgSetup   $trgSetup"
	    #live
	    #./runana 112 $inRootFile $outHistF $nEv_max $rndseed
	    #screen
	    #screenName='tr'$jobID
            #echo "$screenName"
            #screen -S $screenName -L -d -m ./runana 112 $inRootFile $outHistF $nEv_max $rndseed
	    #srun
	    #srun $simHomeDir/runana 112 $inRootFile $outHistF $nEv_max $rndseed
	    if [ "$typeJob" = "sbatch" ]; then
		srun $simHomeDir/runana 112 $inRootFile $outHistF $nEv_max $rndseed $trgSetup
	    elif [ "$typeJob" = "screen" ]; then
		screenName='tr'$jobID
		echo "$screenName"
		screen -S $screenName -L -d -m $simHomeDir/runana 112 $inRootFile $outHistF $nEv_max $rndseed $trgSetup
	    elif [ "$typeJob" = "live" ]; then
		$simHomeDir/runana 112 $inRootFile $outHistF $nEv_max $rndseed $trgSetup
	    else
		printHelp
	    fi
	else
	    printHelp
	fi
    elif [ "$1" = "-test_live_NSB" ]; then
	if [ $# -eq 2 ]; then
	    inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/"
	    nEv_max=$2
	    jobID="0000"
	    inRootFile=$inRootFilePref$jobID"/corsika_"$jobID"ID.root"
	    outHistF="./hist_trgA_corsika_"$jobID"ID_test_live.root"
	    outlogF="./hist_trgA_corsika_"$jobID"ID_test_live.log"
	    rndseed=`date +%N`
	    #rndseed=742827582
	    echo "inRootFile $inRootFile"	
	    echo "outHistF   $outHistF"
	    echo "jobID      $jobID"
	    echo "nEv_max    $nEv_max"
	    echo "rndseed    $rndseed"
	    #live
	    ./runana 112 $inRootFile $outHistF $nEv_max $rndseed | tee $outlogF
	else
	    printHelp
	fi
    elif [ "$1" = "-test_live_NSB_k_dist" ]; then
	if [ $# -eq 1 ]; then
	    inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/"
	    nEv_max=3
	    jobID="0000"
	    inRootFile=$inRootFilePref$jobID"/corsika_"$jobID"ID.root"
	    outHistF="./hist_trgA_corsika_k_dist_"$jobID"ID_test_live.root"
	    outlogF="./hist_trgA_corsika_k_dist_"$jobID"ID_test_live.log"
	    rndseed=`date +%N`
	    #rndseed=742827582
	    echo "inRootFile $inRootFile"	
	    echo "outHistF   $outHistF"
	    echo "jobID      $jobID"
	    echo "nEv_max    $nEv_max"
	    echo "rndseed    $rndseed"
	    #live
	    ./runana 113 $inRootFile $outHistF $nEv_max $rndseed | tee $outlogF
	else
	    printHelp
	fi
    elif [ "$1" = "-test_live" ]; then
	if [ $# -eq 2 ]; then
	    nEv_max=$2
	    inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/"
	    outHistFPref="./"
	    binE=0
	    binTheta=0
	    binDist=0
	    jobID="0000"
	    data_chunk_ID=-999
	    inRootFile=$inRootFilePref$jobID"/corsika_"$jobID"ID.root"
	    mkdir -p $outHistFPref$jobID
	    outHistF=$outHistFPref"/hist_trgA_test_live_proton_corsika_"$jobID"ID.root"
	    outlogF=$outHistFPref"/hist_trgA_test_live_proton_corsika_"$jobID"ID.log"
	    npe_min=20
	    npe_max=10000
	    rndseed=`date +%N`
	    #
	    echo "inRootFile    $inRootFile"
	    echo "outHistF      $outHistF"
	    echo "binE          $binE"
	    echo "binTheta      $binTheta"
	    echo "binDist       $binDist"
	    echo "jobID         $jobID"
	    echo "npe_min       $npe_min"	
	    echo "npe_max       $npe_max"
	    echo "nEv_max       $nEv_max"
	    echo "rndseed       $rndseed"
	    echo "data_chunk_ID $data_chunk_ID"
	    #
	    #live
	    ./runana 111 $inRootFile $outHistF $binE $binTheta $binDist $npe_min $npe_max $nEv_max $rndseed $data_chunk_ID | tee $outlogF
	else
	    printHelp
	fi
    elif [ "$1" = "-test_srun" ]; then
	if [ $# -eq 2 ]; then
	    inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/"
	    nEv_max=$2
	    jobID="0000"
	    inRootFile=$inRootFilePref$jobID"/corsika_"$jobID"ID.root"
	    outHistF="./hist_trgA_corsika_"$jobID"ID_test_live.root"
	    rndseed=`date +%N`
	    echo "inRootFile $inRootFile"	
	    echo "outHistF   $outHistF"
	    echo "jobID      $jobID"
	    echo "nEv_max    $nEv_max"
	    echo "rndseed    $rndseed"
	    #live
	    srun ./runana 112 $inRootFile $outHistF $nEv_max $rndseed
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


