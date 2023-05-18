add_executable(heat_seq main.c solver_seq.c)
target_link_libraries(heat_seq PUBLIC heat_common)

add_executable(heat_ompss2 main.c solver_ompss2.c)
target_link_libraries(heat_ompss2 PRIVATE heat_common
  Nanos6::wrapper)

add_executable(heat_ompss2_residual main.c solver_ompss2_residual.c)
target_link_libraries(heat_ompss2_residual PRIVATE heat_common
  Nanos6::wrapper)
