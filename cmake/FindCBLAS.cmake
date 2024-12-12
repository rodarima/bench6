include(GNUInstallDirs)

find_library(CBLAS_LIBRARY NAMES cblas)
find_path(CBLAS_INCLUDE_DIR cblas.h)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(CBLAS DEFAULT_MSG
  CBLAS_LIBRARY CBLAS_INCLUDE_DIR)

if(NOT CBLAS_FOUND)
  return()
endif()

if(TARGET CBLAS)
  return()
endif()

add_library(CBLAS SHARED IMPORTED)
set_target_properties(CBLAS PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${CBLAS_INCLUDE_DIR}"
  IMPORTED_LOCATION ${CBLAS_LIBRARY})
