include(GNUInstallDirs)

find_library(TAMPI_LIBRARY NAMES tampi-c)
find_path(TAMPI_INCLUDE_DIR TAMPI.h)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Tampi DEFAULT_MSG
  TAMPI_LIBRARY TAMPI_INCLUDE_DIR)

if(NOT TAMPI_FOUND)
  return()
endif()

if(TARGET Tampi::tampi-c)
  return()
endif()

add_library(Tampi::tampi-c SHARED IMPORTED)
set_target_properties(Tampi::tampi-c PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${TAMPI_INCLUDE_DIR}"
  IMPORTED_LOCATION ${TAMPI_LIBRARY})
