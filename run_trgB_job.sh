#!/bin/bash
#SBATCH --job-name crgen%j
#SBATCH --error /srv/beegfs/scratch/users/b/burmistr/pyeventio_example/job_error/crgen_%j.error
#SBATCH --output /srv/beegfs/scratch/users/b/burmistr/pyeventio_example/job_output/output_%j.output
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition public-cpu
#SBATCH --time 0-03:00:00

simHomeDir="./"

source $simHomeDir/setupEnv.sh -d
source $simHomeDir/convert2root.sh --set_unlimited_mem

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d   : default"
    echo " [1]      : particle type (g,gd,e,p,nsb)"
    echo " [2]      : jobID         - [0:199]"
    echo " [3]      : jobType       - [sbatch:screen:live]"
    echo " [0] -c   : recompile"
    echo " [0] -h   : print help"
}

if [ $# -eq 0 ]; then
    printHelp
else
    if [ "$1" = "-d" ]; then
	if [ $# -eq 4 ]; then
	    if [ "$2" = "g" ]; then
		inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/root/"
		outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/trgB/"
		NGBsim=0
		anaConf="./anaTrgB_g.conf"
	    elif [ "$2" = "gd" ]; then
		inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_diffuse_nsb_1x/root/"
		outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_diffuse_nsb_1x/trgB/"
		NGBsim=0
		anaConf="./anaTrgB_gd.conf"
	    elif [ "$2" = "e" ]; then
		inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/electron_nsb_1x/root/"
		outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/electron_nsb_1x/trgB/"
		NGBsim=0
		anaConf="./anaTrgB_e.conf"
	    elif [ "$2" = "p" ]; then
		#inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/"
		#outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/trgB/"
		inRootFilePref="../scratch/simtel_data/proton/root/"
		outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/trgB/"
		NGBsim=0
		anaConf="./anaTrgB_p.conf"
	    elif [ "$2" = "nsb" ]; then
		inRootFilePref="../scratch/simtel_data/proton/root/"
		outHistFPref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgB/"
		NGBsim=1
		anaConf="./anaTrgB_nsb.conf"
	    fi
	    jobID=$3
	    typeJob=$4
	    trgSetup="trg_setup.conf"	    
	    rndseed=`date +%N`
	    #
	    #inRootFile=$inRootFilePref$jobID"/corsika_"$jobID"ID.root"
	    #mkdir -p $outHistFPref$jobID
	    #outHistF=$outHistFPref$jobID"/hist_trgB_corsika_"$jobID"ID.root"
	    inRootFile=$inRootFilePref"/corsika_run"$jobID".root"
	    outHistF=$outHistFPref"/hist_trgB_corsika_run"$jobID".root"
	    npe_min=0
	    npe_max=1000000
	    nEvSim_max=-1
	    #nEv_max=15000
	    #nEv_max=-1
	    nEv_max=100
            NGB_rate_in_MHz="265.0"
	    #NGB_rate_in_MHz="0.1"
	    fadc_electronic_noise_RMS="3.8082498"
	    #
	    #
	    echo "inRootFile                $inRootFile"
	    echo "outHistF                  $outHistF"
	    echo "NGB_rate_in_MHz           $NGB_rate_in_MHz"
	    echo "fadc_electronic_noise_RMS $fadc_electronic_noise_RMS"
	    echo "npe_min                   $npe_min"	
	    echo "npe_max                   $npe_max"
	    echo "nEvSim_max                $nEvSim_max"
	    echo "nEv_max                   $nEv_max"
	    echo "rndseed                   $rndseed"
	    echo "NGBsim                    $NGBsim"
	    echo "trgSetup                  $trgSetup"
	    echo "anaConf                   $anaConf"
	    #
	    #
	    if [ "$typeJob" = "sbatch" ]; then
		srun $simHomeDir/runana 222 $inRootFile $outHistF $NGB_rate_in_MHz $fadc_electronic_noise_RMS $npe_min $npe_max $nEvSim_max $nEv_max $rndseed $NGBsim $trgSetup $anaConf
	    elif [ "$typeJob" = "screen" ]; then
		screenName='tr'$jobID
		echo "$screenName"
		screen -S $screenName -L -d -m $simHomeDir/runana 222 $inRootFile $outHistF $NGB_rate_in_MHz $fadc_electronic_noise_RMS $npe_min $npe_max $nEvSim_max $nEv_max $rndseed $NGBsim $trgSetup  $anaConf
	    elif [ "$typeJob" = "live" ]; then
		#echo "$typeJob"
		$simHomeDir/runana 222 $inRootFile $outHistF $NGB_rate_in_MHz $fadc_electronic_noise_RMS $npe_min $npe_max $nEvSim_max $nEv_max $rndseed $NGBsim $trgSetup $anaConf
	    else
		printHelp
	    fi
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


