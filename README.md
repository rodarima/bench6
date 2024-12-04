# Bench6 - A benchmark suite

This benchmark suite is a monorepo that contains a set of micro-benchmarks
for OmpSs-2 and several mini-apps.

Each mini-app is available in different variants that include multiple
programming models. The CMake build system is used to find the available
dependencies and enable them accordingly. All the binaries are prefixed with
`b6_` so they can be easily found with auto-complete.

    b6_heat_nanos6 -s 2048 -t 10 -b 64
    b6_nbody_omp -f -b 512
    ...

To install all the variants, you will need:

- Clang with support for OpenMP and OpenMP-V as well as OmpSs-2
- OpenMP
- OpenMP-V
- Nanos6
- NODES
- nOS-V
- MPI
- TAMPI
- GPI-2
- Tagaspi
- BLAS
- LAPACK

A Nix package is available that includes all the variants. Use `nix develop` to
build the benchmarks and enter a shell where you can run them.

The bigotes tool can be used to run any mini-app several times, so that enough
executions are done to have a representative measurement. Example:

    bigotes b6_heat_nanos6 -s 2048 -t 10 -b 64
