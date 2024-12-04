#!/bin/sh

set -e
set -x

#B=bigotes

$B b6_heat_nanos6 -s 2048 -t 10 -b 64
$B b6_nbody_omp -f -b 512
$B b6_nbody_ompv -f -b 512
$B b6_miniamr_mpi_omp
$B b6_nqueens_nanos6
$B b6_strassen_nanos6
$B b6_cholesky_nanos6
