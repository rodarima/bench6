include(GNUInstallDirs)

find_library(NOSV_LIBRARY NAMES nosv PATHS "$ENV{NOSV_HOME}" "$ENV{NOSV_HOME}/lib")
find_path(NOSV_INCLUDE_DIR nosv.h PATHS "$ENV{NOSV_HOME}" "$ENV{NOSV_HOME}/include")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Nosv DEFAULT_MSG
  NOSV_LIBRARY NOSV_INCLUDE_DIR)

if(NOT NOSV_FOUND)
  return()
endif()

if(TARGET Nosv::nosv)
  return()
endif()

add_library(Nosv::nosv SHARED IMPORTED)
set_target_properties(Nosv::nosv PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${NOSV_INCLUDE_DIR}"
  IMPORTED_LOCATION ${NOSV_LIBRARY})
