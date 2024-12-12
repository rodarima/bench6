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
