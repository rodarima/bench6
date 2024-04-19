get_property(BENCH6_LIST GLOBAL PROPERTY bench6_list)
configure_file(config.in.h config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(bigote runner.c)
target_link_libraries(bigote PRIVATE m bench6_lib)
install(TARGETS bigote RUNTIME DESTINATION bin)
