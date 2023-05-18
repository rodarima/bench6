add_library(bench6_common STATIC common.c)

macro(mk_bench6 NAME SOURCE)
  add_executable(${NAME} ${SOURCE})
  target_link_libraries(${NAME} PRIVATE bench6_common Nanos6::wrapper)
  install(TARGETS ${NAME} RUNTIME DESTINATION bin)
endmacro()

mk_bench6(bench6_register_deps register_deps.c)
mk_bench6(bench6_sched_add sched_add.c)
mk_bench6(bench6_sched_get sched_get.c)
mk_bench6(bench6_readywave readywave.c)
