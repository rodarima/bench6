set(COMMON lulesh-viz.cc lulesh-util.cc)

if(NANOS6_FOUND)
  mk_bench(b6_lulesh_nanos6_task)
  target_sources(b6_lulesh_nanos6_task PRIVATE ${COMMON}
    lulesh-init.cc lulesh-comm.cc lulesh_t.cc)
  target_compile_definitions(b6_lulesh_nanos6_task
    PRIVATE "USE_MPI=0" "NANOS6_RUNTIME")
  target_link_libraries(b6_lulesh_nanos6_task
    PRIVATE Nanos6::nanos6)
endif()
