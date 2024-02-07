#!/bin/sh

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d    : default"
    echo " [1]       : particle type (g,gd,e,p)"
    echo " [0] -h    : print help"
}

counter=0

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	if [ $# -eq 2 ]; then
	    particletype=$2
	    if [ "$particletype" = "g" ]; then
                inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/root/"
                nDir=10
		n_multiply=25
		log_out="multiply_data_gamma_on_nsb_1x.log"
            elif [ "$particletype" = "gd" ]; then
                inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_diffuse_nsb_1x/root/"
                nDir=20
		n_multiply=75
		log_out="multiply_data_gamma_diffuse_nsb_1x.log"
            elif [ "$particletype" = "e" ]; then
                inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/electron_nsb_1x/root/"
                nDir=20
		n_multiply=10
		log_out="multiply_data_electron_nsb_1x.log"
            elif [ "$particletype" = "p" ]; then
                inRootFilePref="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/"
                nDir=50
		n_multiply=30
		log_out="multiply_data_proton_nsb_1x.log"
	    fi
            for i in `seq 0 $(echo "$nDir-1" | bc )` ; do
		folderID=`printf %04d $i`
		file_to_multiply=$inRootFilePref/$folderID/"corsika_"$folderID"ID.root"
		for j in `seq 0 $(echo "$n_multiply-1" | bc )` ; do
		    new_folderID=`printf %04d $(echo "$folderID+($j+1)*$nDir" | bc)`
		    new_file="corsika_"$new_folderID"ID.root"
		    cmd="mkdir $inRootFilePref/$new_folderID"
		    echo "$cmd"
		    $cmd
		    cmd="cp $file_to_multiply $inRootFilePref/$new_folderID/$new_file"
		    echo "$cmd"
		    $cmd | tee $log_out
		done
	    done
	else
	    printHelp
	fi
    elif [ "$1" = "-h" ]; then
	printHelp
    else
        printHelp
    fi
fi

#espeak "I have done"
