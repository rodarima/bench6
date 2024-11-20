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

find_package(OpenMP)
if (OPENMP_FOUND)
  mk_bench(b6_nbody_omp_plain)
  target_sources(b6_nbody_omp_plain PRIVATE ${PLAIN_SOURCES} src/plain/solver_omp.c)
  target_link_libraries(b6_nbody_omp_plain PRIVATE OpenMP::OpenMP_C)
endif()
