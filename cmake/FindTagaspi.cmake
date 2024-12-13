include(GNUInstallDirs)

find_library(TAGASPI_LIBRARY NAMES tagaspi PATHS "$ENV{TAGASPI_HOME}" "$ENV{TAGASPI_HOME}/lib")
find_path(TAGASPI_INCLUDE_DIR TAGASPI.h PATHS "$ENV{TAGASPI_HOME}" "$ENV{TAGASPI_HOME}/include")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Tagaspi DEFAULT_MSG
  TAGASPI_LIBRARY TAGASPI_INCLUDE_DIR)

if(NOT TAGASPI_FOUND)
  return()
endif()

if(TARGET Tagaspi::tagaspi)
  return()
endif()

add_library(Tagaspi::tagaspi SHARED IMPORTED)
set_target_properties(Tagaspi::tagaspi PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${TAGASPI_INCLUDE_DIR}"
  IMPORTED_LOCATION ${TAGASPI_LIBRARY})
