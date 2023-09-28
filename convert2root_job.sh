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
    echo " [3]       : run_i"
    echo " [0] -h    : print help"
}

function creat_csv2toot_list {
    #corsika_run1.header.pkl_0.csv
    #corsika_run1.pe_info.pkl_0.csv
    csv_dir=$1
    csv_dir_int=$2
    n_run_header=`ls $csv_dir/corsika_run*.header.pkl_0.csv | wc -w`
    n_run_pe=`ls $csv_dir/corsika_run*.pe_info.pkl_0.csv | wc -w`
    list="$csv_dir/list"
    rm -rf $list
    run_i_start=$(echo "$csv_dir_int*100 + 1" | bc)
    run_i_stop=$(echo "$csv_dir_int*100 + $n_run_header" | bc)
    echo " creat_csv2toot_list     --> $csv_dir"
    echo " n_run_header            --> $n_run_header"
    echo " n_run_pe                --> $n_run_pe"
    echo " list                    --> $list"
    echo " run_i_start             --> $run_i_start"
    echo " run_i_stop              --> $run_i_stop"
    #
    for run_i in `seq $run_i_start $run_i_stop`
    do
	n_files_header=`ls $csv_dir/corsika_run$run_i.header.pkl_*.csv | wc -w`
	n_files_pe=`ls $csv_dir/corsika_run$run_i.pe_info.pkl_*.csv | wc -w`
	echo " run: $run_i n_files_header   --> $n_files_header"
	echo "         n_run_pe         --> $n_files_pe"
	for pkl_i in `seq 0 $(echo "$n_run_header-1" | bc )`
	do
	    header_csv_file="corsika_run$run_i.header.pkl_$pkl_i.csv"
	    pe_info_csv_file="corsika_run$run_i.pe_info.pkl_$pkl_i.csv"	    
	    if [ -f "$csv_dir/$header_csv_file" ]; then
		if [ -f "$csv_dir/$pe_info_csv_file" ]; then
		    #echo "$header_csv_file $pe_info_csv_file" | tee -a $csv_dir/list
		    echo "$header_csv_file $pe_info_csv_file" >> $list
		fi
		#echo "$header_csv_file" | tee -a $list
		#echo "$header_csv_file" >> $list
	    fi
	done
    done
}

function load_modules {
    module load GCC/7.3.0-2.30 GCCcore/7.3.0 OpenMPI/3.1.1 ROOT/6.14.06-Python-2.7.15 Geant4/10.5 CMake/3.11.4
}
    
if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
        if [ $# -eq 4 ]; then
	    load_modules
	    ulimit -s unlimited
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
	    folderID=$3
	    i=$4
	    csv_dir="/srv/beegfs/scratch/users/b/burmistr/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/csv/$folderID/"
	    outputRootFilePath="/srv/beegfs/scratch/users/b/burmistr/mono-lst-sipm-pmma-3ns-v1_triggerless/$particlein/root/$folderID/"
	    mkdir -p $outputRootFilePath
	    #creat_csv2toot_list $csv_dir $i
	    srun /home/users/b/burmistr/pyeventio_example/convert2root 2 $csv_dir "$csv_dir/list" $outputRootFilePath/"corsika_"$folderID"ID.root"
	else
            printHelp
        fi
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
