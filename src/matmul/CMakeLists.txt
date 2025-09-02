include_directories(src)

if(NANOS6_FOUND AND CBLAS_FOUND)
  mk_bench(b6_matmul_nanos6)
  target_sources(b6_matmul_nanos6 PRIVATE src/smp/matmul_ompss2.c)
  target_link_libraries(b6_matmul_nanos6 PRIVATE CBLAS Nanos6::nanos6)
endif()

if(NODES_FOUND AND CBLAS_FOUND)
  mk_bench(b6_matmul_nodes)
  target_sources(b6_matmul_nodes PRIVATE src/smp/matmul_ompss2.c)
  target_link_libraries(b6_matmul_nodes PRIVATE CBLAS Nodes::nodes)
endif()

if(OPENMP_FOUND AND CBLAS_FOUND)
  mk_bench(b6_matmul_omp_iter)
  target_sources(b6_matmul_omp_iter PRIVATE src/smp/matmul_openmp_iter.c)
  target_link_libraries(b6_matmul_omp_iter PRIVATE CBLAS OpenMP::OpenMP_C)
endif()

if(OPENMPV_FOUND AND CBLAS_FOUND)
  mk_bench(b6_matmul_ompv_iter)
  target_sources(b6_matmul_ompv_iter PRIVATE src/smp/matmul_openmp_iter.c)
  target_link_libraries(b6_matmul_ompv_iter PRIVATE CBLAS OpenMPV)
endif()

if(MPI_FOUND AND TAMPI_FOUND AND CBLAS_FOUND)
  set(MPI_SRC src/common/matmul.c src/mpi/main.c)
  if(NANOS6_FOUND)
    mk_bench(b6_matmul_tampi_nanos6)
    target_sources(b6_matmul_tampi_nanos6 PRIVATE ${MPI_SRC} src/mpi/01.matmul_ompss2_tampi.c)
    target_link_libraries(b6_matmul_tampi_nanos6
      PRIVATE MPI::MPI_C Tampi::tampi-c CBLAS Nanos6::nanos6)

    mk_bench(b6_matmul_itampi_nanos6)
    target_sources(b6_matmul_itampi_nanos6 PRIVATE ${MPI_SRC} src/mpi/01.matmul_ompss2_tampi.c)
    target_link_libraries(b6_matmul_itampi_nanos6
      PRIVATE MPI::MPI_C Tampi::tampi-c CBLAS Nanos6::nanos6)
  endif()

  if(NODES_FOUND)
    mk_bench(b6_matmul_tampi_nodes)
    target_sources(b6_matmul_tampi_nodes PRIVATE ${MPI_SRC} src/mpi/01.matmul_ompss2_tampi.c)
    target_link_libraries(b6_matmul_tampi_nodes
      PRIVATE MPI::MPI_C Tampi::tampi-c CBLAS Nodes::nodes)

    mk_bench(b6_matmul_itampi_nodes)
    target_sources(b6_matmul_itampi_nodes PRIVATE ${MPI_SRC} src/mpi/02.matmul_ompss2_itampi.c)
    target_link_libraries(b6_matmul_itampi_nodes
      PRIVATE MPI::MPI_C Tampi::tampi-c CBLAS Nodes::nodes)
  endif()
endif()

# Macro to set up MKL, TAMPI, and MPI configuration for targets
macro(mk_bench_mpi_tampi_mkl target)
  target_compile_definitions(${target} PRIVATE "USE_MKL")
  target_link_libraries(${target} PRIVATE MPI::MPI_C Tampi::tampi-c MKL::MKL)
endmacro()

if(MPI_FOUND AND TAMPI_FOUND AND TARGET MKL::MKL)
  set(MPI_SRC src/common/matmul.c src/mpi/main.c)
  if(NANOS6_FOUND)
    mk_bench(b6_matmul_tampi_mkl_nanos6)
    target_sources(b6_matmul_tampi_mkl_nanos6 PRIVATE ${MPI_SRC} src/mpi/01.matmul_ompss2_tampi.c)
    mk_bench_mpi_tampi_mkl(b6_matmul_tampi_mkl_nanos6)
    target_link_libraries(b6_matmul_tampi_mkl_nanos6 PRIVATE Nanos6::nanos6)

    mk_bench(b6_matmul_itampi_mkl_nanos6)
    target_sources(b6_matmul_itampi_mkl_nanos6 PRIVATE ${MPI_SRC} src/mpi/01.matmul_ompss2_tampi.c)
    mk_bench_mpi_tampi_mkl(b6_matmul_itampi_mkl_nanos6)
    target_link_libraries(b6_matmul_itampi_mkl_nanos6 PRIVATE Nanos6::nanos6)
  endif()

  if(NODES_FOUND)
    mk_bench(b6_matmul_tampi_mkl_nodes)
    target_sources(b6_matmul_tampi_mkl_nodes PRIVATE ${MPI_SRC} src/mpi/01.matmul_ompss2_tampi.c)
    mk_bench_mpi_tampi_mkl(b6_matmul_tampi_mkl_nodes)
    target_link_libraries(b6_matmul_tampi_mkl_nodes PRIVATE Nodes::nodes)

    mk_bench(b6_matmul_itampi_mkl_nodes)
    target_sources(b6_matmul_itampi_mkl_nodes PRIVATE ${MPI_SRC} src/mpi/02.matmul_ompss2_itampi.c)
    mk_bench_mpi_tampi_mkl(b6_matmul_itampi_mkl_nodes)
    target_link_libraries(b6_matmul_itampi_mkl_nodes PRIVATE Nodes::nodes)
  endif()

  if (NOSV_FOUND)
    mk_bench(b6_matmul_tampi_mkl_nosv)
    target_compile_definitions(b6_matmul_tampi_mkl_nosv PRIVATE USE_NOSV_ONLY)
    target_sources(b6_matmul_tampi_mkl_nosv PRIVATE ${MPI_SRC} src/mpi/03.matmul_nosv_tampi.c)
    mk_bench_mpi_tampi_mkl(b6_matmul_tampi_mkl_nosv)
    target_link_libraries(b6_matmul_tampi_mkl_nosv PRIVATE Nosv::nosv)
  endif()
endif()

