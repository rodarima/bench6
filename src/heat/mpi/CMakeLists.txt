if(NOT MPI_FOUND)
  return()
endif()

macro(mk_heat_mpi NAME SOURCE)
  mk_bench(${NAME})
  target_sources(${NAME} PRIVATE ${SOURCE})
  target_link_libraries(${NAME} PRIVATE heat_mpi_common)
endmacro()

macro(mk_heat_mpi_nanos6 NAME SOURCE)
  mk_heat_mpi(${NAME} ${SOURCE})
  target_compile_options(${NAME} PRIVATE "-fompss-2=libnanos6")
  target_link_options(${NAME} PRIVATE "-fompss-2=libnanos6")
endmacro()

macro(mk_heat_mpi_nodes NAME SOURCE)
  mk_heat_mpi(${NAME} ${SOURCE})
  target_link_libraries(${NAME} PRIVATE Nodes::wrapper)
endmacro()

macro(mk_heat_tampi_nanos6 NAME SOURCE)
  mk_heat_mpi_nanos6(${NAME} ${SOURCE})
  target_link_libraries(${NAME} PRIVATE Tampi::tampi-c)
endmacro()


# -------------------------------------------------------------------

add_library(heat_mpi_common STATIC main.c utils.c)
target_link_libraries(heat_mpi_common PUBLIC heat_common MPI::MPI_C)

mk_heat_mpi(b6_heat_mpi solver_mpi.c)
mk_heat_mpi(b6_heat_mpi_nbuffer solver_mpi_nbuffer.c)

if(NANOS6_FOUND)
  mk_heat_mpi_nanos6(b6_heat_mpi_nanos6_forkjoin solver_mpi_ompss2_forkjoin.c)
  mk_heat_mpi_nanos6(b6_heat_mpi_nanos6_tasks solver_mpi_ompss2_tasks.c)
endif()

if(NODES_FOUND)
  mk_heat_mpi_nodes(b6_heat_mpi_nodes_forkjoin solver_mpi_ompss2_forkjoin.c)
  mk_heat_mpi_nodes(b6_heat_mpi_nodes_tasks solver_mpi_ompss2_tasks.c)
endif()

if(NANOS6_FOUND AND TAMPI_FOUND)
  mk_heat_tampi_nanos6(b6_heat_itampi_nanos6_tasks solver_itampi_ompss2_tasks.c)
endif()
