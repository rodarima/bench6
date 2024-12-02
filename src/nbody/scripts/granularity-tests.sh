#!/bin/bash

# Check if Makefile is in the default folder
cd ${0%/*} # Never assume the directory and always change it to where the script is located

makefile_path=".."
if [ ! -f "${makefile_path}/Makefile" ]; then
    read -p 'Makefile not found by default, please specify in which folder it is found (relative path): ' makefile_path
fi

function usage {
    echo "Usage: ./granularity-tests.sh <PARAMETERS> [-bs blocksize] [-bo BIGO] [OPTIONAL_PARAMETERS]\n"
    echo "Parameters:\n"
    echo "  -p PARTICLES: use PARTICLES as the total number of particles\n"
    echo "  -t TIMESTEPS: use TIMESTEPS as the number of timesteps\n\n"
    echo "Optional parameters:\n"
    echo "  -bo BO - use BO as the BIGO (defaults to N2)\n"
    echo "  -bs BS - use BS as blocksize (defaults to 2048)\n"
    echo "  -f, - forcefully generate particles by avoiding creating files\n"
    echo "  -c, - check the correctness of the result (disabled by default)\n"
    echo "  -C, - do not check the correctness of the result\n"
    echo "  -o, - save the computed particles to the default output file (disabled by default)\n"
    echo "  -O, - do not save the computed particles to the default output file\n"
    echo "  -h, - display this help and exit\n\n"
    exit 1
}


#####################
# COLLECT ARGUMENTS #
#####################

argc=$#
argv=($@)

for (( j=0; j<argc; j++ )); do
    opt=${argv[$j]}

    if [ "$opt" == "-bs" ]; then
        ((j++))
        blocksize="${argv[$j]}"
    elif [ "$opt" == "-bo" ]; then
        ((j++))
        bigo="${argv[$j]}"
    elif [ "$opt" == "-p" ] || [ "$opt" == "-t" ]; then
        ((j++))
        parameters="$parameters $opt ${argv[$j]}"
    elif [ "$opt" == "-f" ] || [ "$opt" == "-c" ] || [ "$opt" == "-C" ] || [ "$opt" == "-o" ] || [ "$opt" == "-O" ]; then
        parameters="$parameters $opt"
    elif [ "$opt" == "-h" ]; then
        usage
    fi
done


###########################
# CALL THE CORRECT BINARY #
###########################

# Compile
cd $makefile_path

if [ "$blocksize" != "" ] && [ "$bigo" != "" ]; then
    BS=$blocksize BIGO=$bigo make
elif [ "$blocksize" != "" ]; then
    BS=$blocksize make
    bigo="N2"
elif [ "$bigo" != "" ]; then
    BIGO=$bigo make
    blocksize="2048"
else
    blocksize="2048"
    bigo="N2"
    make
fi

# Execute

# By default use the ompss version
#output=$(./nbody_mpi.omp.${blocksize}.exe $parameters)
#output=$(./nbody_mpi_ompss.${bigo}.${blocksize}.exe $parameters)
#output=$(./nbody_omp_plain.${bigo}.exe $parameters)
 output=$(./nbody_ompss.${bigo}.${blocksize}.exe $parameters)
#output=$(./nbody_seq.${bigo}.${blocksize}.exe $parameters)
#output=$(./nbody_seq_plain.${bigo}.exe $parameters)

echo $output


