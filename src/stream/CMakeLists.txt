# all: stream_gcc stream_ompss2 stream_intel
# 
# CFLAGS=-Wall -mcmodel=medium
# 
# stream_gcc:
# 	gcc -O3 -fopenmp -ffast-math $(CFLAGS) kernel_openmp.c main_openmp.c -o stream_gcc
# 
# stream_intel:
# 	icc -O3 -qopenmp $(CFLAGS) kernel_openmp.c main_openmp.c -o stream_intel
# 
# stream_ompss2:
# 	mcc -O3 --ompss-2 --Wl,-ffast-math $(CFLAGS) kernel_ompss2.c main_ompss2.c -o stream_ompss2
# 
# clean:
# 	rm -f stream_gcc stream_ompss2 stream_intel

if(NANOS6_FOUND)
  mk_bench(b6_stream_nanos6)
  target_sources(b6_stream_nanos6 PRIVATE kernel_ompss2.c main.c)
  target_link_libraries(b6_stream_nanos6 PRIVATE Nanos6::nanos6)
endif()

if(OPENMP_FOUND)
  mk_bench(b6_stream_omp)
  target_sources(b6_stream_omp PRIVATE kernel_openmp.c main.c)
  target_link_libraries(b6_stream_omp PRIVATE OpenMP::OpenMP_C)
endif()
