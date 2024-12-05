# Streaming

## Introduction
The Streaming is a communication-intensive benchmark that taskifies both computations and
communications. The participating computing nodes process multiple large chunks of data,
and each node has a different function that must apply to each chunk. There are no data
dependencies between elements of a chunk when applying a function. From the first node, a
chunk moves across all the rest of the nodes, one by one, applying each function to the
chunk's elements. That process is repeated for each chunk, one after the other, following
a pipeline of chunk computations across the nodes. The chunks are computed and transferred
by blocks. The blocks of a chunk can be processed by a node in parallel. Processes avoid
intra-node communication and always transfer data to other nodes. Processes have a receive
and a send buffer with sufficient space to receive/send all the blocks of a single chunk
simultaneously.

## Software requirements

The Streaming application requires the following software:

  * [OmpSs-2](https://github.com/bsc-pm/ompss-2-releases) installation with the Nanos6 runtime and LLVM/Clang compiler supporting OmpSs-2.
  * MPI installation with support for MPI_THREAD_MULTIPLE (`mpicc` or `mpiicc`).
  * [TAMPI](https://github.com/bsc-pm/tampi) installation. The `TAMPI_HOME` envar must be defined to its installation directory. *Optional*.
  * [TAGASPI](https://github.com/bsc-pm/tagaspi) installation. The `TAGASPI_HOME` envar must be defined to its installation directory. The TAGASPI variants
    require a special [GPI-2](https://pm.bsc.es/gitlab/interoperability/extern/GPI-2) installation (branch `lowlevel`) with support for the lowlevel operations
    extension. The `GASPI_EXT_HOME` envar must be defined to its installation directory. *Optional*.

## Available versions and building instructions

The Streaming application has several versions which are compiled in different
binaries, by executing the `make` command. They are:

  * **01.streaming_tampi_ompss2_tasks.bin**: Parallel version using MPI + OmpSs-2 tasks + Blocking TAMPI.
  * **02.streaming_itampi_ompss2_tasks.bin**: Parallel version using MPI + OmpSs-2 tasks + Non-blocking TAMPI.
  * **03.streaming_tagaspi_ompss2_tasks.bin**: Parallel version using GASPI + OmpSs-2 tasks + TAGASPI. Acks for one-sided are waited by extra tasks.
  * **04.streaming_tagaspi_ompss2_tasks_onready.bin**: Parallel version using GASPI + OmpSs-2 tasks + TAGASPI. Acks for one-sided are waited with onready clause.
  * **05.streaming_mpi.bin**: Parallel version using MPI and non-blocking primitives.


  The simplest way to compile this package is:

  1. Stay in the Streaming root directory to recursively build all the versions.
     The TAMPI versions are compiled only if the environment variable
     `TAMPI_HOME` is set to the Task-Aware MPI (TAMPI) library installation
     directory. Similarly, the TAGASPI versions are compiled if `TAGASPI_HOME`
     and `GASPI_EXT_HOME` are defined to the Task-Aware GASPI (TAGASPI) library
     and the extended GPI-2 installation directories.

     Note that TAGASPI versions need to be linked to a GASPI implementation
     (`GASPI_EXT_HOME`) supporting the lowlevel operations extension, which
     is a new API used by the TAGASPI library to privide the GASPI+OmpSs-2
     interoperability features.

  2. Type `make` to compile the selected benchmark's version(s). The default
     Makefile is for the Clang compiler (Makefile.clang).

## Execution instructions

The binaries accept several options. The mandatory options are the global size,
the block size, and the number of timesteps to perform. The global size is determined
by the `-s SIZE` option, which indicates the global array size across all process. The
chunk size is the size of array that each process will have and where all the
operations will be applied. The chunk size is determined by dividing the global size
(SIZE) by the number of processes.

Additionally, the chunks are processed and transfered by blocks. The block size is
specified by the `-b BSIZE` option. The block size is the granularity of computations,
communications, and tasks (if applicable). The `-t TSTEPS` option specifies the number
of timesteps. The timesteps indicate how many times each process will apply the
operations on the blocks of the chunk.

Finally, the `-o OFFSET` option is useful when the program is run by more than one
process per node and you want to guarantee that communication is always inter-node.
This option can be used to avoid intra-node communication between consecutive processes
that are placed within the same node. The OFFSET specifies how to compute the neighboring
processes for each process. A process with rank R will receive chunk data from the
process with rank R-OFFSET and will send chunk data to the process with rank R+OFFSET.
By default, the OFFSET is one. The help `-h` option shows other optional parameters.

An example of execution with MPICH could be:

```
$ mpiexec.hydra -n 4 -bind-to hwthread:16 ./01.streaming_tampi_ompss2_tasks.bin -s 8192 -b 128 -t 150
```

Or with the SLURM launcher:

```
$ srun --cpu_bind=cores -n 4 -c 16 ./01.streaming_tampi_ompss2_tasks.bin -s 8192 -b 128 -t 150
```

in which the application will perform 150 timesteps in 4 MPI processes with 16 hardware
threads (or cores) in each process that will be used by the OmpSs-2 runtime. The size
of the global array is 8192. Thus, the chunk size for each process is 8192 / 4 processes =
2048 elements. The block size is 128, so each process will divide its chunk by 16 blocks.
