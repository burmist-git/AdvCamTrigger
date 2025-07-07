#!/bin/bash -l

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d                : default"
    echo " [2]                   : Nnodes (1-...)"
    echo " [0] --copyForAllTypes : make separate folders for all particle types"
    echo " [0] -c                : clean"
    echo " [0] -h                : print help"
}

nCPU_per_node=72
nCPU_idle=1
nJOB_per_node=$(echo "$nCPU_per_node - $nCPU_idle" | bc -l)
greasyJobDir="./csv_pyeventio_stereo_greasy_job/"
outGreasySbatch_sh="./run_csv_pyeventio_stereo_greasy_sbatch.sh"
mkdir -p $greasyJobDir

if [ $# -eq 0 ]; then
    printHelp
else
    if [ "$1" = "-d" ]; then
        if [ $# -eq 2 ]; then
	    #
	    #particletype="gamma"
	    #particletype="gamma_diffuse"
	    #particletype="electron"
	    #particletype="proton"
	    #particletype="NSB268MHz"
	    #
	    particletype="proton"
	    #inFilePref="/scratch/snx3000/lburmist/simtel_data/proton_st/data/"
	    inFilePref="../scratch/simtel_data/gamma/data/"
	    #outFilePref="/scratch/snx3000/lburmist/simtel_data/proton_st/csv/"
	    outFilePref="../scratch/simtel_data/gamma/csv/"
	    #
	    Nnodes=$2
	    echo "Nnodes        $Nnodes"
	    echo "nCPU_per_node $nCPU_per_node"
	    echo "nCPU_idle     $nCPU_idle"
	    echo "nJOB_per_node $nJOB_per_node"
	    echo "particletype  $particletype"
	    echo "inFilePref    $inFilePref"
	    echo "outFilePref   $outFilePref"
	    echo "greasyJobDir       $greasyJobDir"
	    echo "outGreasySbatch_sh $outGreasySbatch_sh"
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
		    mkdir -p $outFilePref
		    echo "  fileID      = $fileID"
		    echo "  inFilePref  = $inFilePref"
		    echo "  outFilePref = $outFilePref"		    
		    #
		    #/scratch/snx3000/lburmist/simtel_data/$particletype/data/corsika_run1.simtel.gz
		    in_simtel_file="$inFilePref/corsika_run$fileID.simtel.gz"
		    #
		    #if [ -f "$in_simtel_file" ]; then
			#
			out_header_file="$outFilePref/corsika_run$fileID.header.csv"
			out_pe_inf_file_LST1="$outFilePref/corsika_run$fileID.pe_info_LST1.csv"
			out_pe_inf_file_LST2="$outFilePref/corsika_run$fileID.pe_info_LST2.csv"
			out_pe_inf_file_LST3="$outFilePref/corsika_run$fileID.pe_info_LST3.csv"
			out_pe_inf_file_LST4="$outFilePref/corsika_run$fileID.pe_info_LST4.csv"			
			#
			rm -rf $out_header_file
			rm -rf $out_pe_inf_file_LST1
			rm -rf $out_pe_inf_file_LST2
			rm -rf $out_pe_inf_file_LST3
			rm -rf $out_pe_inf_file_LST4
			#			
			cmd="python csv_pyeventio_stereo.py $in_simtel_file $out_header_file $out_pe_inf_file_LST1 $out_pe_inf_file_LST2 $out_pe_inf_file_LST3 $out_pe_inf_file_LST4"
			echo "$cmd" >> $outJOBfileList
		    #fi
		    #
		    ((fileID=fileID+1))
		done
	    done
	else
	    printHelp   
	fi      
    elif [ "$1" = "-c" ]; then
	rm -rf $greasyJobDir/*
	rm -rf $outGreasySbatch_sh
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
