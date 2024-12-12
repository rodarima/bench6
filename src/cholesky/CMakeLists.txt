if(CBLAS_FOUND AND LAPACKE_FOUND)
  if(NANOS6_FOUND)
    mk_bench(b6_cholesky_nanos6)
    target_sources(b6_cholesky_nanos6 PRIVATE cholesky.c)
    target_link_libraries(b6_cholesky_nanos6
      PRIVATE Nanos6::nanos6 CBLAS LAPACKE)
  endif()

  if(NODES_FOUND)
    mk_bench(b6_cholesky_nodes)
    target_sources(b6_cholesky_nodes PRIVATE cholesky.c)
    target_link_libraries(b6_cholesky_nodes
      PRIVATE Nodes::nodes CBLAS LAPACKE)
  endif()
endif()
