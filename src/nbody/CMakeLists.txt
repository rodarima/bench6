include_directories(src)

set(PLAIN_SOURCES
  src/common/common.c
  src/plain/utils.c
  src/plain/main.c)

set(SMP_SOURCES
  src/common/common.c
  src/blocking/common/common_utils.c
  src/blocking/smp/utils.c
  src/blocking/smp/main.c)

mk_bench(b6_nbody_seq_plain)
target_sources(b6_nbody_seq_plain PRIVATE ${PLAIN_SOURCES} src/plain/solver_seq.c)
target_link_libraries(b6_nbody_seq_plain PRIVATE m)

mk_bench(b6_nbody_seq)
target_sources(b6_nbody_seq PRIVATE ${SMP_SOURCES} src/blocking/smp/solver_seq.c)
target_link_libraries(b6_nbody_seq PRIVATE m)

if (OPENMP_FOUND)
  mk_bench(b6_nbody_omp_plain)
  target_sources(b6_nbody_omp_plain PRIVATE ${PLAIN_SOURCES} src/plain/solver_omp.c)
  target_link_libraries(b6_nbody_omp_plain PRIVATE m OpenMP::OpenMP_C)

  mk_bench(b6_nbody_omp)
  target_sources(b6_nbody_omp PRIVATE ${SMP_SOURCES} src/blocking/smp/solver_omp.c)
  target_link_libraries(b6_nbody_omp PRIVATE m OpenMP::OpenMP_C)
endif()

if (OPENMPV_FOUND)
  mk_bench(b6_nbody_ompv)
  target_sources(b6_nbody_ompv PRIVATE ${SMP_SOURCES} src/blocking/smp/solver_omp.c)
  target_link_libraries(b6_nbody_omp PRIVATE m OpenMPV)
endif()
