add_subdirectory(core)

include_directories(core)

if(NANOS6_FOUND)
  mk_bench(b6_taskbench_nanos6)
  target_sources(b6_taskbench_nanos6 PRIVATE ompss/main.cc)
  target_link_libraries(b6_taskbench_nanos6 PRIVATE taskbench_core Nanos6::nanos6)
  target_compile_definitions(b6_taskbench_nanos6 PRIVATE "OMPSS2")
endif()

if(NODES_FOUND)
  mk_bench(b6_taskbench_nodes)
  target_sources(b6_taskbench_nodes PRIVATE ompss/main.cc)
  target_link_libraries(b6_taskbench_nodes PRIVATE taskbench_core Nodes::nodes)
  target_compile_definitions(b6_taskbench_nodes PRIVATE "OMPSS2")
endif()

if(OPENMP_FOUND)
  mk_bench(b6_taskbench_omp)
  target_sources(b6_taskbench_omp PRIVATE openmp/main.cc)
  target_link_libraries(b6_taskbench_omp PRIVATE taskbench_core OpenMP::OpenMP_C)
endif()

if(OPENMPV_FOUND)
  mk_bench(b6_taskbench_ompv)
  target_sources(b6_taskbench_ompv PRIVATE openmp/main.cc)
  target_link_libraries(b6_taskbench_ompv PRIVATE taskbench_core OpenMPV)
endif()
