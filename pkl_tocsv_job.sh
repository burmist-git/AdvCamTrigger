#!/bin/sh
#SBATCH --job-name easc%j
#SBATCH --error /srv/beegfs/scratch/users/b/burmistr/pyeventio_example/job_error/crgen_%j.error
#SBATCH --output /srv/beegfs/scratch/users/b/burmistr/pyeventio_example/job_output/output_%j.output
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition public-cpu
#SBATCH --time 24:00:00

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d    : single job"
    echo " [1]       : part_type (g/gd/e/p)"
    echo " [2]       : jobID"
    echo " [0] -h    : print help"
}

function conv_to_csv {
    #
    #corsika_run1.header.pkl_0
    #corsika_run1.pe_info.pkl_0
    #../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/pkl/0000/
    #
    particle=$1
    particlein=$2
    folderID=$3
    pklpath="/srv/beegfs/scratch/users/b/burmistr/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/pkl/$folderID/"
    csvpath="/srv/beegfs/scratch/users/b/burmistr/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/$folderID/"
    mkdir -p $csvpath
    #
    echo "particle   $particle"
    echo "particlein $particlein"
    echo "folderID   $folderID"
    echo "pklpath    $pklpath"
    echo "csvpath    $csvpath"
    #
    for filename_header_pkl in `ls -lrt $pklpath/corsika_run*.header.pkl_* | awk {'print $9'}`
    do
        filename_header_csv=$csvpath`basename $filename_header_pkl`'.csv'
	#ls $filename_header_pkl
        srun python pkl_tocsv.py header $filename_header_pkl $filename_header_csv
    done
    #
    for filename_pe_pkl in `ls -lrt $pklpath/corsika_run*.pe_info.pkl_* | awk {'print $9'}`
    do
        filename_pe_csv=$csvpath`basename $filename_pe_pkl`'.csv'
	#ls $filename_pe_pkl
        srun python pkl_tocsv.py pe $filename_pe_pkl $filename_pe_csv
    done
}

if [ $# -eq 0 ] 
then
    printHelp
else
    if [ "$1" = "-d" ]; then
        if [ $# -eq 3 ]; then
            if [ "$2" = "g" ]; then
                particle="gamma"
                particlein="gamma_on_nsb_1x"
            elif [ "$2" = "gd" ]; then
                particle="gamma_diffuse"
                particlein="gamma_diffuse_nsb_1x"
            elif [ "$2" = "e" ]; then
                particle="electron"
                particlein="electron_nsb_1x"
            elif [ "$2" = "p" ]; then
                particle="proton"
                particlein="proton_nsb_1x"
            else
                printHelp
            fi
	    jobID=$3
	    conv_to_csv $particle $particlein $jobID
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
