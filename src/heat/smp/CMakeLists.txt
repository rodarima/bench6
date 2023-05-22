macro(mk_heat_smp NAME SOURCE)
  mk_bench(${NAME})
  target_sources(${NAME} PRIVATE ${SOURCE})
  target_link_libraries(${NAME} PRIVATE heat_smp_common)
endmacro()

macro(mk_heat_nanos6 NAME SOURCE)
  mk_heat_smp(${NAME} ${SOURCE})
  target_link_libraries(${NAME} PRIVATE Nanos6::wrapper)
endmacro()

macro(mk_heat_nodes NAME SOURCE)
  mk_heat_smp(${NAME} ${SOURCE})
  target_link_libraries(${NAME} PRIVATE Nodes::wrapper)
endmacro()

# -------------------------------------------------------------------

add_library(heat_smp_common STATIC main.c)
target_link_libraries(heat_smp_common PUBLIC heat_common)

mk_heat_smp(b6_heat_seq solver_seq.c)

if(NANOS6_FOUND)
  mk_heat_nanos6(b6_heat_nanos6 solver_ompss2.c)
  mk_heat_nanos6(b6_heat_nanos6_residual solver_ompss2_residual.c)
endif()

if(NODES_FOUND)
  mk_heat_nodes(b6_heat_nodes solver_ompss2.c)
  mk_heat_nodes(b6_heat_nodes_residual solver_ompss2_residual.c)
endif()
