if(NANOS6_FOUND AND MPI_FOUND)
  if(TAMPI_FOUND)
      mk_bench(b6_tsunampi_tampi)
      target_sources(b6_tsunampi_tampi PRIVATE main.c tampi.c)
      target_link_libraries(b6_tsunampi_tampi
        PRIVATE m MPI::MPI_C Tampi::tampi-c Nanos6::nanos6)
  endif()
  if(GASPI_FOUND AND TAGASPI_FOUND)
      mk_bench(b6_tsunampi_tagaspi)
      target_sources(b6_tsunampi_tagaspi PRIVATE main.c tagaspi.c)
      target_link_libraries(b6_tsunampi_tagaspi
        PRIVATE m MPI::MPI_C Gaspi::gaspi Tagaspi::tagaspi Nanos6::nanos6)
  endif()
endif()
