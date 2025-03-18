
# CPPFLAGS = -I. -DOMP_EXTENSIONS
# CFLAGS   = -fopenmp -Ofast -march=native -ffast-math
# LDFLAGS  = -fopenmp -Ofast -march=native -ffast-math


# MPI is needed for all targets
if(NOT MPI_FOUND)
  return()
endif()

add_compile_options("-DOMP_EXTENSIONS")

set(SOURCES block.c check_sum.c comm_block.c comm.c comm_parent.c comm_refine.c
  comm_util.c def.c driver.c init.c move.c pack.c plot.c profile.c rcb.c
  refine.c stencil.c util.c)

add_library(miniamr STATIC ${SOURCES})

if (OPENMP_FOUND)
  mk_bench(b6_miniamr_mpi_omp)
  target_sources(b6_miniamr_mpi_omp PRIVATE main.c)
  target_link_libraries(b6_miniamr_mpi_omp PRIVATE miniamr m MPI::MPI_C OpenMP::OpenMP_C)
endif()

if (OPENMPV_FOUND)
  mk_bench(b6_miniamr_mpi_ompv)
  target_sources(b6_miniamr_mpi_ompv PRIVATE main.c)
  target_link_libraries(b6_miniamr_mpi_ompv PRIVATE miniamr m MPI::MPI_C)
  target_compile_options(b6_miniamr_mpi_ompv PRIVATE "-fopenmp=libompv")
  target_link_options(b6_miniamr_mpi_ompv PRIVATE "-fopenmp=libompv")
endif()
