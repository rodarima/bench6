#!/bin/sh

set -e
set -x

bigotes b6_heat_nanos6 -s 2048 -t 10 -b 64
bigotes b6_nbody_omp -f -b 512
bigotes b6_nbody_ompv -f -b 512
