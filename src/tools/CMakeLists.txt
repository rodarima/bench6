get_property(BENCH6_LIST GLOBAL PROPERTY bench6_list)
configure_file(config.in.h config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(bench6_runner runner.c)
target_link_libraries(bench6_runner PRIVATE m bench6_lib)
install(TARGETS bench6_runner RUNTIME DESTINATION bin)
