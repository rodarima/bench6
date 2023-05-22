add_library(heat_common STATIC common/misc.c common/kernel.c)
target_include_directories(heat_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(heat_common PUBLIC m)

add_subdirectory(smp)
add_subdirectory(mpi)
