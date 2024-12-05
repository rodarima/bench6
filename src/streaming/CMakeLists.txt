include_directories(src)


if(MPI_FOUND)
  set(MPI_SRC src/common/misc.c src/mpi/main.c src/common/kernel.c)

  mk_bench(b6_streaming_mpi)
  target_sources(b6_streaming_mpi PRIVATE ${MPI_SRC} src/mpi/05.solver_mpi.c)
  target_link_libraries(b6_streaming_mpi PRIVATE m MPI::MPI_C)

  if(NANOS6_FOUND)
    if(TAMPI_FOUND)
      mk_bench(b6_streaming_tampi_nanos6)
      target_sources(b6_streaming_tampi_nanos6
        PRIVATE ${MPI_SRC} src/mpi/01.solver_tampi_ompss2_tasks.c)
      target_link_libraries(b6_streaming_tampi_nanos6
        PRIVATE m MPI::MPI_C Tampi::tampi-c Nanos6::nanos6)
  
      mk_bench(b6_streaming_itampi_nanos6)
      target_sources(b6_streaming_itampi_nanos6
        PRIVATE ${MPI_SRC} src/mpi/02.solver_itampi_ompss2_tasks.c)
      target_link_libraries(b6_streaming_itampi_nanos6
        PRIVATE m MPI::MPI_C Tampi::tampi-c Nanos6::nanos6)
    endif()
  
    if(GASPI_FOUND AND TAGASPI_FOUND)
      set(GASPI_SRC src/common/misc.c src/gaspi/utils.c src/gaspi/main.c
        src/common/kernel.c)
  
      mk_bench(b6_streaming_tagaspi_nanos6)
      target_sources(b6_streaming_tagaspi_nanos6
        PRIVATE ${GASPI_SRC} src/gaspi/03.solver_tagaspi_ompss2_tasks.c)
      target_link_libraries(b6_streaming_tagaspi_nanos6
        PRIVATE m MPI::MPI_C Gaspi::gaspi Tagaspi::tagaspi Nanos6::nanos6)
  
      mk_bench(b6_streaming_tagaspi_nanos6_onready)
      target_sources(b6_streaming_tagaspi_nanos6_onready
        PRIVATE ${GASPI_SRC} src/gaspi/04.solver_tagaspi_ompss2_tasks_onready.c)
      target_link_libraries(b6_streaming_tagaspi_nanos6_onready
        PRIVATE m MPI::MPI_C Gaspi::gaspi Tagaspi::tagaspi Nanos6::nanos6)
    endif()
  endif(NANOS6_FOUND)
endif(MPI_FOUND)
