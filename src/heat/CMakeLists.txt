# --- Common to all ---

add_library(heat_kernel STATIC kernel.c)
target_include_directories(heat_kernel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(heat_kernel PRIVATE 
	-Rpass-analysis=loop-vectorize
	-ffast-math)

add_library(heat_common STATIC misc.c)
target_include_directories(heat_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(heat_common PUBLIC m heat_kernel)

# --- SMP based ---

add_library(heat_smp_common STATIC main_smp.c)
target_link_libraries(heat_smp_common PUBLIC heat_common)

macro(mk_heat_smp NAME SOURCE)
  mk_bench(${NAME})
  target_sources(${NAME} PRIVATE ${SOURCE})
  target_link_libraries(${NAME} PRIVATE heat_smp_common)
endmacro()

# No requisites
mk_heat_smp(b6_heat_seq solver_seq.c)

if(NANOS6_FOUND)
  macro(mk_heat_nanos6 NAME SOURCE)
    mk_heat_smp(${NAME} ${SOURCE})
    target_link_libraries(${NAME} PRIVATE Nanos6::nanos6)
  endmacro()

  mk_heat_nanos6(b6_heat_nanos6           solver_ompss2.c)
  mk_heat_nanos6(b6_heat_nanos6_residual  solver_ompss2_residual.c)
endif()

if(NODES_FOUND)
  macro(mk_heat_nodes NAME SOURCE)
    mk_heat_smp(${NAME} ${SOURCE})
    target_link_libraries(${NAME} PRIVATE Nodes::nodes)
  endmacro()

  mk_heat_nodes(b6_heat_nodes             solver_ompss2.c)
  mk_heat_nodes(b6_heat_nodes_residual    solver_ompss2_residual.c)
endif()

if (OPENMP_FOUND)
  mk_heat_smp(b6_heat_omp solver_openmp.c)
  target_link_libraries(b6_heat_omp PRIVATE OpenMP::OpenMP_C)
endif()

if (OPENMPV_FOUND)
  mk_heat_smp(b6_heat_ompv solver_openmp.c)
  target_link_libraries(b6_heat_ompv PRIVATE OpenMPV)
endif()

# --- MPI based ---

if(MPI_FOUND)
  macro(mk_heat_mpi NAME SOURCE)
    mk_bench(${NAME})
    target_sources(${NAME} PRIVATE ${SOURCE})
    target_link_libraries(${NAME} PRIVATE heat_mpi_common)
  endmacro()

  add_library(heat_mpi_common STATIC main_mpi.c utils_mpi.c)
  target_link_libraries(heat_mpi_common PUBLIC heat_common MPI::MPI_C)

  mk_heat_mpi(b6_heat_mpi solver_mpi.c)
  mk_heat_mpi(b6_heat_mpi_nbuffer solver_mpi_nbuffer.c)

  if(NANOS6_FOUND)
    macro(mk_heat_mpi_nanos6 NAME SOURCE)
      mk_heat_mpi(${NAME} ${SOURCE})
      target_link_libraries(${NAME} PRIVATE Nanos6::nanos6)
    endmacro()
    mk_heat_mpi_nanos6(b6_heat_mpi_nanos6_forkjoin  solver_mpi_ompss2_forkjoin.c)
    mk_heat_mpi_nanos6(b6_heat_mpi_nanos6_tasks     solver_mpi_ompss2_tasks.c)
  endif()

  if(NODES_FOUND)
    macro(mk_heat_mpi_nodes NAME SOURCE)
      mk_heat_mpi(${NAME} ${SOURCE})
      target_link_libraries(${NAME} PRIVATE Nodes::nodes)
    endmacro()
    mk_heat_mpi_nodes(b6_heat_mpi_nodes_forkjoin solver_mpi_ompss2_forkjoin.c)
    mk_heat_mpi_nodes(b6_heat_mpi_nodes_tasks solver_mpi_ompss2_tasks.c)
  endif()

  if(OPENMP_FOUND)
    mk_heat_mpi(b6_heat_mpi_omp solver_mpi_openmp.c)
    target_link_libraries(b6_heat_mpi_omp PRIVATE OpenMP::OpenMP_C)
  endif()

  if(OPENMPV_FOUND)
    mk_heat_mpi(b6_heat_mpi_ompv solver_mpi_openmp.c)
    target_link_libraries(b6_heat_mpi_ompv PRIVATE OpenMPV)

  endif()

  # --- TAMPI based (requires MPI) ---

  if (TAMPI_FOUND)
    macro(mk_heat_tampi NAME SOURCE)
      mk_heat_mpi(${NAME} ${SOURCE})
      target_link_libraries(${NAME} PRIVATE Tampi::tampi-c)
    endmacro()

    if(NANOS6_FOUND)
      mk_heat_tampi(b6_heat_itampi_nanos6_tasks solver_itampi_ompss2_tasks.c)
      target_link_libraries(b6_heat_itampi_nanos6_tasks PRIVATE Nanos6::nanos6)
    endif()

    if(NODES_FOUND)
      mk_heat_tampi(b6_heat_itampi_nodes_tasks solver_itampi_ompss2_tasks.c)
      target_link_libraries(b6_heat_itampi_nodes_tasks PRIVATE Nodes::nodes)
    endif()

    if(OPENMPV_FOUND)
      mk_heat_tampi(b6_heat_tampi_ompv solver_tampi_openmp.c)
      target_link_libraries(b6_heat_tampi_ompv PRIVATE OpenMPV)
    endif(OPENMPV_FOUND)

  endif(TAMPI_FOUND)
endif(MPI_FOUND)

# --- GASPI ---

# TODO
