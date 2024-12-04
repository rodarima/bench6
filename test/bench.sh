#!/bin/sh

set -e
set -x

bigotes b6_heat_nanos6 -s 2048 -t 10 -b 64
bigotes b6_nbody_omp -f -b 512
bigotes b6_nbody_ompv -f -b 512
bigotes b6_miniamr_mpi_omp
bigotes b6_nqueens_nanos6
bigotes b6_strassen_nanos6
