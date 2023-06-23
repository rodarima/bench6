add_library(heat_kernel STATIC common/kernel.c)
target_include_directories(heat_kernel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(heat_kernel PRIVATE 
	-Rpass-analysis=loop-vectorize
	-ffast-math)

add_library(heat_common STATIC common/misc.c)
target_include_directories(heat_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(heat_common PUBLIC m heat_kernel)

add_subdirectory(smp)
add_subdirectory(mpi)
