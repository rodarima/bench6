
# CPPFLAGS = -I. -DOMP_EXTENSIONS
# CFLAGS   = -fopenmp -Ofast -march=native -ffast-math
# LDFLAGS  = -fopenmp -Ofast -march=native -ffast-math


# MPI and TAMPI is needed for all targets
if(MPI_FOUND AND TAMPI_FOUND)

  set(SOURCES block.c check_sum.c comm_block.c comm.c comm_parent.c
    comm_refine.c comm_util.c def.c driver.c init.c move.c pack.c plot.c
    profile.c rcb.c refine.c stencil.c kernel.c util.c main.c)

  # We need to rebuild all sources for Nanos6 and Nodes, as the compiler will
  # create different code to interface with the runtimes.

  if(NANOS6_FOUND)
    mk_bench(b6_miniamr_tampi_nanos6)
    target_sources(b6_miniamr_tampi_nanos6 PRIVATE ${SOURCES})
    target_link_libraries(b6_miniamr_tampi_nanos6
      PRIVATE m MPI::MPI_C Tampi::tampi-c Nanos6::nanos6)
    target_compile_definitions(b6_miniamr_tampi_nanos6 PRIVATE "USE_NANOS6")
  endif()

  if(NODES_FOUND)
    mk_bench(b6_miniamr_tampi_nodes)
    target_sources(b6_miniamr_tampi_nodes PRIVATE ${SOURCES})
    target_link_libraries(b6_miniamr_tampi_nodes
      PRIVATE m MPI::MPI_C Tampi::tampi-c Nodes::nodes)
    target_compile_definitions(b6_miniamr_tampi_nodes PRIVATE "USE_NODES")
  endif()
endif()
