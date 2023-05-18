if(NOT MPI_FOUND)
  return()
endif()

add_library(heat_mpi_common STATIC main.c utils.c)
target_link_libraries(heat_mpi_common PUBLIC heat_common MPI::MPI_C)

macro(mk_heat_mpi NAME SOURCE)
  add_executable(${NAME} ${SOURCE})
  target_link_libraries(${NAME} PRIVATE heat_mpi_common)
  install(TARGETS ${NAME} RUNTIME DESTINATION bin)
endmacro()

mk_heat_mpi(heat_mpi solver_mpi.c)
mk_heat_mpi(heat_mpi_nbuffer solver_mpi_nbuffer.c)

if(NANOS6_FOUND)
  macro(mk_heat_mpi_nanos6 NAME SOURCE)
    mk_heat_mpi(${NAME} ${SOURCE})
    target_link_libraries(${NAME} PRIVATE Nanos6::wrapper)
  endmacro()

  mk_heat_mpi_nanos6(heat_mpi_nanos6_forkjoin solver_mpi_ompss2_forkjoin.c)
  mk_heat_mpi_nanos6(heat_mpi_nanos6_tasks solver_mpi_ompss2_tasks.c)
endif()

if(NODES_FOUND)
  macro(mk_heat_mpi_nodes NAME SOURCE)
    mk_heat_mpi(${NAME} ${SOURCE})
    target_link_libraries(${NAME} PRIVATE Nodes::wrapper)
  endmacro()

  mk_heat_mpi_nodes(heat_mpi_nodes_forkjoin solver_mpi_ompss2_forkjoin.c)
  mk_heat_mpi_nodes(heat_mpi_nodes_tasks solver_mpi_ompss2_tasks.c)
endif()
