set(SOURCES box.c comm.c comm_send.c comm_recv.c comm_pack.c comm_unpack.c
  comm_wait.c force.c input.c integrate.c ljs.c neigh.c packbuf.c packbuf_mpi.c
  ref.c safe.c thermo.c)

if(MPI_FOUND AND TAMPI_FOUND)
  if(NANOS6_FOUND)
    mk_bench(b6_minimd_tampi_nanos6)
    target_sources(b6_minimd_tampi_nanos6 PRIVATE ${SOURCES})
    target_link_libraries(b6_minimd_tampi_nanos6
      PRIVATE m MPI::MPI_C Tampi::tampi-c Nanos6::nanos6)
  endif()
  if(NODES_FOUND)
    mk_bench(b6_minimd_tampi_nodes)
    target_sources(b6_minimd_tampi_nodes PRIVATE ${SOURCES})
    target_link_libraries(b6_minimd_tampi_nodes
      PRIVATE m MPI::MPI_C Tampi::tampi-c Nodes::nodes)
  endif()
endif()
