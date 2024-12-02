#!/bin/bash

run_verbose() {
   echo "=== Running: ${@} ==="
   echo ""
   eval "${@}"
}

# Assuming MPICH and 2 nodes of 64-core processor with single NUMA
# Leveraging 4 x 4 x 4 blocks per node
# Testing rank/cores configurations:
#  - 1  ranks/node: 2 x 1 x 1 ranks, 4 x 4 x 4 blocks, 8 x 4 x 4 = 128 total blocks
#  - 2  ranks/node: 2 x 2 x 1 ranks, 4 x 2 x 4 blocks, 8 x 4 x 4 = 128 total blocks
#  - 4  ranks/node: 2 x 2 x 2 ranks, 4 x 2 x 2 blocks, 8 x 4 x 2 = 128 total blocks
#  - 8  ranks/node: 4 x 2 x 2 ranks, 2 x 2 x 2 blocks, 8 x 4 x 4 = 128 total blocks
#  - 16 ranks/node: 4 x 4 x 2 ranks, 2 x 1 x 2 blocks, 8 x 4 x 4 = 128 total blocks
# Matching with the MPI-only configuration:
#  - 64 ranks/node: 8 x 4 x 4 ranks, 1 x 1 x 1 blocks, 8 x 4 x 4 = 128 total blocks

# Number of configurations to test
nconfs=5

# System configuration
nodes=2
coresxnode=64

# Block size
bs=12

# Number of ranks per dimension
npx=(2 2 2 4 4)
npy=(1 2 2 2 4)
npz=(1 1 2 2 2)

# Number of blocks per rank per dimension
nbx=(4 4 4 2 2)
nby=(4 2 2 2 1)
nbz=(4 4 2 2 2)

# Maximum number of blocks (approximate)
maxblocks=(55000 30000 18000 10000 7000)

# Extra parameters
extra_params="--send_faces --separate_buffers --max_comm_tasks 8"

# Input objects
obj1="--object 2 0 0.25 0.25 0.25 0.016 0.0 0.0 0.23 0.23 0.23 0.0 0.0 0.0"
obj2="--object 2 0 0.25 0.75 0.75 0.016 0.0 0.0 0.23 0.23 0.23 0.0 0.0 0.0"
obj3="--object 2 0 0.75 0.25 0.75 -0.016 0.0 0.0 0.23 0.23 0.23 0.0 0.0 0.0"
obj4="--object 2 0 0.75 0.75 0.25 -0.016 0.0 0.0 0.23 0.23 0.23 0.0 0.0 0.0"
objects="--num_objects 4 $obj1 $obj2 $obj3 $obj4"

# Setting OmpSs-2 variables
export NANOS6=optimized
export NANOS6_PRIORITY=0
export NANOS6_DEPENDENCIES=discrete

for (( c=0; c<$nconfs; c++ )); do
   # Computing execution parameters
   ranks=$((${npx[$c]}*${npy[$c]}*${npz[$c]}))
   coresxrank=$((($nodes*$coresxnode)/$ranks))
   ranksxnode=$(($ranks/$nodes))

   # Setting OpenMP variables
   export OMP_NUM_THREADS=$coresxrank

   echo "Running configuration $c with $ranks ranks and $coresxrank cores/rank"
   for rep in {1..1}; do
      run_verbose mpiexec.hydra -n $ranks -bind-to core:$coresxrank -ppn $ranksxnode ./miniAMR.x --npx ${npx[$c]} --npy ${npy[$c]} --npz ${npz[$c]} --nx $bs --ny $bs --nz $bs --init_x ${nbx[$c]} --init_y ${nby[$c]} --init_z ${nbz[$c]} --max_blocks ${maxblocks[$c]} --num_vars 40 --stages_per_ts 40 --checksum_freq 10 --stencil 7 --num_tsteps 14 --refine_freq 5 --num_refine 4 $extra_params $objects
   done
done
