set(SOURCES 
  src/mpi_interop_ompss/HPCCG.cpp
  src/common/ddot.cpp
  src/common/read_HPC_row.cpp
  src/common/compute_residual.cpp
  src/common/dump_matlab_matrix.cpp
  src/common/main.cpp
  src/common/YAML_Doc.cpp
  src/common/HPC_sparsemv.cpp
  src/common/mytimer.cpp
  src/common/HPC_Sparse_Matrix.cpp
  src/common/waxpby.cpp
  src/common/generate_matrix.cpp
  src/common/YAML_Element.cpp)

include_directories(include)
add_compile_definitions(WALL)

if(OPENMP_FOUND)
  mk_bench(b6_hpccg_omp_task)
  target_sources(b6_hpccg_omp_task PRIVATE ${SOURCES})
  target_link_libraries(b6_hpccg_omp_task PRIVATE m OpenMP::OpenMP_C)
endif()

if(OPENMPV_FOUND)
  mk_bench(b6_hpccg_ompv_task)
  target_sources(b6_hpccg_ompv_task PRIVATE ${SOURCES})
  target_link_libraries(b6_hpccg_ompv_task PRIVATE m OpenMPV)
endif()
