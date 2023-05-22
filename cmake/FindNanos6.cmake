include(GNUInstallDirs)

if(DEFINED ENV{NANOS6_HOME})
  set(NANOS6_HOME "$ENV{NANOS6_HOME}")
else()
  message(STATUS "NANOS6_HOME not set, refusing to search")
endif()

find_library(NANOS6_LIBRARY NAMES nanos6 PATHS "${NANOS6_HOME}/lib" NO_DEFAULT_PATH)
find_path(NANOS6_INCLUDE_DIR nanos6.h PATHS "${NANOS6_HOME}/include" NO_DEFAULT_PATH)
find_file(NANOS6_WRAPPER NAMES nanos6-main-wrapper.o PATHS "${NANOS6_HOME}/lib" NO_DEFAULT_PATH)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Nanos6 DEFAULT_MSG
  NANOS6_LIBRARY NANOS6_INCLUDE_DIR NANOS6_WRAPPER)

if(NOT NANOS6_FOUND)
  return()
endif()

if(TARGET Nanos6::nanos6)
  return()
endif()

add_library(Nanos6::nanos6 SHARED IMPORTED)
set_target_properties(Nanos6::nanos6 PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${NANOS6_INCLUDE_DIR}"
  IMPORTED_LOCATION ${NANOS6_LIBRARY})

add_library(Nanos6::wrapper STATIC IMPORTED)
set_target_properties(Nanos6::wrapper PROPERTIES
  IMPORTED_LOCATION ${NANOS6_WRAPPER})
target_compile_options(Nanos6::wrapper INTERFACE "-fompss-2")
target_link_libraries(Nanos6::wrapper INTERFACE Nanos6::nanos6)
