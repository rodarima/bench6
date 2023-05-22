add_library(bench6_lib STATIC bench6.c common.c)
target_include_directories(bench6_lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
