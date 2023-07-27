#!/bin/sh

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d                 : setup environment"
    echo " [0] -p                 : purge modules"
    echo " [0] -list              : list modules"
    echo " [0] -h                 : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	module load GCC/7.3.0-2.30 GCCcore/7.3.0 OpenMPI/3.1.1 ROOT/6.14.06-Python-2.7.15 Geant4/10.5 CMake/3.11.4
	#module load GCC/9.3.0 GSL/2.6
	#module load GCC/9.3.0 GSL/2.6 OpenMPI/4.0.3 SciPy-bundle/2020.03-Python-3.8.2
    elif [ "$1" = "-p" ]; then
	module purge
    elif [ "$1" = "-list" ]; then
	module list
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
