add_subdirectory(core)

include_directories(core)

if(NANOS6_FOUND)
  mk_bench(b6_taskbench_nanos6)
  target_sources(b6_taskbench_nanos6 PRIVATE ompss/main.cc)
  target_link_libraries(b6_taskbench_nanos6 PRIVATE taskbench_core Nanos6::nanos6)
  target_compile_definitions(b6_taskbench_nanos6 PRIVATE "OMPSS2")
endif()
