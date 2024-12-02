include(GNUInstallDirs)

find_library(GASPI_LIBRARY NAMES GPI2)
find_path(GASPI_INCLUDE_DIR GASPI.h)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Gaspi DEFAULT_MSG
  GASPI_LIBRARY GASPI_INCLUDE_DIR)

if(NOT GASPI_FOUND)
  return()
endif()

if(TARGET Gaspi::Gaspi)
  return()
endif()

add_library(Gaspi::gaspi SHARED IMPORTED)
set_target_properties(Gaspi::gaspi PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${GASPI_INCLUDE_DIR}"
  IMPORTED_LOCATION ${GASPI_LIBRARY})
