#!/bin/bash

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] --createdata : create data sample"
    echo " [0] -h           : print help"
}

function gen_data {
    #
    #particle="gamma"
    #particlein="gamma_on_nsb_1x"
    #for f in `ls ../scratch/data_nagaia/data/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/output/log/corsika_run*.log`
    #do
    #	grep "@\*" $f > $particle/`basename $f`
    #done
    #tar -czvf $particle.tar.gz $particle
    #
    particle=$1
    particlein=$2
    echo "$particle"
    echo "$particlein"
    for f in `ls ../scratch/data_nagaia/data/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/output/log/corsika_run*.log`
    do
	echo $f
	grep "@\*" $f > $particle/`basename $f`
    done
    tar -czvf $particle.tar.gz $particle
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "--createdata" ]; then
	#gen_data gamma gamma_on_nsb_1x
	gen_data gamma_diffuse gamma_diffuse_nsb_1x
	gen_data electron electron_nsb_1x
	gen_data proton proton_nsb_1x
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
