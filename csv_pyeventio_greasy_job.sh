#!/bin/bash -l

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d                : default"
    echo " [1]                   : particle type (g,gd,e,p,NSB)"
    echo " [2]                   : Nnodes (1-...)"
    echo " [0] --copyForAllTypes : make separate folders for all particle types"
    echo " [0] -c                : clean"
    echo " [0] -h                : print help"
}

nCPU_per_node=72
nCPU_idle=1
nJOB_per_node=$(echo "$nCPU_per_node - $nCPU_idle" | bc -l)
greasyJobDir="./greasy_job/"
outGreasySbatch_sh="./run_greasy_sbatch.sh"

function copyForAllTypes {
    #
    mkdir -p greasy_job_gamma/$greasyJobDir
    mkdir -p greasy_job_gamma_diffuse/$greasyJobDir
    mkdir -p greasy_job_electron/$greasyJobDir
    mkdir -p greasy_job_proton/$greasyJobDir
    mkdir -p greasy_job_NSB/$greasyJobDir
    #
    cp csv_pyeventio.py greasy_job_gamma/.
    cp csv_pyeventio_greasy_job.sh greasy_job_gamma/.
    #
    cp csv_pyeventio.py greasy_job_gamma_diffuse/.
    cp csv_pyeventio_greasy_job.sh greasy_job_gamma_diffuse/.
    #
    cp csv_pyeventio.py greasy_job_electron/.
    cp csv_pyeventio_greasy_job.sh greasy_job_electron/.
    #
    cp csv_pyeventio.py greasy_job_proton/.
    cp csv_pyeventio_greasy_job.sh greasy_job_proton/.
    #
    cp csv_pyeventio.py greasy_job_NSB/.
    cp csv_pyeventio_greasy_job.sh greasy_job_NSB/.
}

if [ $# -eq 0 ]; then
    printHelp
else
    if [ "$1" = "-d" ]; then
        if [ $# -eq 3 ]; then
            if [ "$2" = "g" ]; then
		particletype="gamma"
	    elif [ "$2" = "gd" ]; then
		particletype="gamma_diffuse"
            elif [ "$2" = "e" ]; then
		particletype="electron"
            elif [ "$2" = "p" ]; then
		particletype="proton"
            elif [ "$2" = "NSB" ]; then
		#particletype="NSB386MHz"
		particletype="NSB268MHz"
            fi
	    #
	    Nnodes=$3
	    echo "Nnodes        $Nnodes"
	    echo "nCPU_per_node $nCPU_per_node"
	    echo "nCPU_idle     $nCPU_idle"
	    echo "nJOB_per_node $nJOB_per_node"
	    #
	    fileID=1
	    #
	    rm -rf $outGreasySbatch_sh
	    #
	    for nodesID in $(seq 1 $Nnodes)
	    do
		#
		echo "nodesID = $nodesID"
		echo "#!/bin/sh" >> $outGreasySbatch_sh
		#
		outJOBfile="$greasyJobDir/node_$nodesID.job"
		outJOBfileList="$greasyJobDir/node_$nodesID.joblist"
		#
		echo "sbatch $outJOBfile" >> $outGreasySbatch_sh
		#
		rm -rf $outJOBfile
		rm -rf $outJOBfileList
		echo "#!/bin/bash -l" >> $outJOBfile
		echo "#SBATCH --job-name=simtel" >> $outJOBfile
		echo "#SBATCH --output=/scratch/snx3000/lburmist/simtel_data/job_outlog/simtel.%j.out" >> $outJOBfile
		echo "#SBATCH --error=/scratch/snx3000/lburmist/simtel_data/job_error/simtel.%j.err" >> $outJOBfile
		echo "#SBATCH --account=cta03" >> $outJOBfile
		echo "#SBATCH --time=24:00:00" >> $outJOBfile
		echo "#SBATCH --nodes=1" >> $outJOBfile
		echo "#SBATCH --cpus-per-task=1" >> $outJOBfile
		echo "#SBATCH --partition=normal" >> $outJOBfile
		echo "#SBATCH --constraint=mc" >> $outJOBfile
		echo " " >> $outJOBfile
		echo "module load daint-mc" >> $outJOBfile
		echo "module load GREASY" >> $outJOBfile
		echo " " >> $outJOBfile
		echo "greasy $outJOBfileList" >> $outJOBfile
		#
		for jobIT in $(seq 1 $nJOB_per_node)
		do
		    #
		    inFilePref="/scratch/snx3000/lburmist/simtel_data/$particletype/data/"
		    outFilePref="/scratch/snx3000/lburmist/simtel_data/$particletype/csv/"
		    mkdir -p $outFilePref
		    echo "  fileID      = $fileID"
		    echo "  inFilePref  = $inFilePref"
		    echo "  outFilePref = $outFilePref"		    
		    #
		    #/scratch/snx3000/lburmist/simtel_data/$particletype/data/corsika_run1.simtel.gz
		    in_simtel_file="$inFilePref/corsika_run$fileID.simtel.gz"
		    #
		    if [ -f "$in_simtel_file" ]; then
			#
			out_header_file="$outFilePref/corsika_run$fileID.header.csv"
			out_pe_inf_file="$outFilePref/corsika_run$fileID.pe_info.csv"
			#
			rm -rf $out_header_file
			rm -rf $out_pe_inf_file
			#
			cmd="python csv_pyeventio.py $in_simtel_file $out_header_file $out_pe_inf_file"
			echo "$cmd" >> $outJOBfileList
		    fi
		    #
		    ((fileID=fileID+1))
		done
	    done
	else
	    printHelp   
	fi      
    elif [ "$1" = "--copyForAllTypes" ]; then
	copyForAllTypes	
    elif [ "$1" = "-c" ]; then
	rm -rf $greasyJobDir/*
	rm -rf $outGreasySbatch_sh
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
