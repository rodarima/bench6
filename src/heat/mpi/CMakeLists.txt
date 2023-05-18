if(NOT MPI_FOUND)
  return()
endif()

macro(mk_heat_mpi NAME SOURCE)
  add_executable(${NAME} ${SOURCE} main.c utils.c)
  target_link_libraries(${NAME} PRIVATE heat_common MPI::MPI_C)
  install(TARGETS ${NAME} RUNTIME DESTINATION bin)
endmacro()

mk_heat_mpi(heat_mpi solver_mpi.c)
mk_heat_mpi(heat_mpi_nbuffer solver_mpi_nbuffer.c)

macro(mk_heat_mpi_oss NAME SOURCE)
  mk_heat_mpi(${NAME} ${SOURCE})
  target_link_libraries(${NAME} PRIVATE Nanos6::wrapper)
endmacro()

mk_heat_mpi_oss(heat_mpi_ompss2_forkjoin solver_mpi_ompss2_forkjoin.c)
mk_heat_mpi_oss(heat_mpi_ompss2_tasks solver_mpi_ompss2_tasks.c)
