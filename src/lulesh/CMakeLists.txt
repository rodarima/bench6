set(COMMON lulesh-viz.cc lulesh-util.cc)


#LULESH_EXEC = lulesh_omp lulesh_t lulesh_tf lulesh_oss_tl

#MPIIMCXX_CXXFLAGS = -g -O3 --ompss-2 --openmp-compatibility -I. -Wall -Wno-comment -Wno-maybe-uninitialized -DUSE_MPI=$(USE_MPI) -DOMPSS -k #--pass-through
#MPIIMCXX_LDFLAGS = -g -O3 --ompss-2 --openmp-compatibility -lstdc++

if(MPI_FOUND AND TAMPI_FOUND)
  add_compile_definitions("USE_MPI=1" "USE_TAMPI=1")

  if(NANOS6_FOUND)
    mk_bench(b6_lulesh_tampi_nanos6_task)
    target_sources(b6_lulesh_tampi_nanos6_task PRIVATE ${COMMON}
      lulesh-init.cc lulesh-comm.cc lulesh_t.cc)
    target_link_libraries(b6_lulesh_tampi_nanos6_task
      PRIVATE MPI::MPI_C Tampi::tampi-c Nanos6::nanos6)
  endif()
endif()
