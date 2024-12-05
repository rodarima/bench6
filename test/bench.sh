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
$B b6_multisaxpy_nanos6
$B b6_streaming_tampi_nanos6 -s 8096 -b 1024 -t 2000
$B b6_matmul_itampi_nanos6 1024 2 512 0

# Output not compatible with bigotes
b6_tsunampi_tampi -r 10
