include(GNUInstallDirs)

find_library(TAMPI_LIBRARY NAMES tampi-c HINTS "$ENV{TAMPI_HOME}" "$ENV{TAMPI_HOME}/lib")
find_path(TAMPI_INCLUDE_DIR TAMPI.h HINTS "$ENV{TAMPI_HOME}" "$ENV{TAMPI_HOME}/include")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Tampi DEFAULT_MSG
  TAMPI_LIBRARY TAMPI_INCLUDE_DIR)

if(NOT TAMPI_FOUND)
  return()
endif()

if(TARGET Tampi::tampi-c)
  return()
endif()

# Parse TAMPI version from the TAMPI_Decl.h file
file(STRINGS "${TAMPI_INCLUDE_DIR}/TAMPI_Decl.h" version-major-line
  REGEX "#define TAMPI_VERSION_MAJOR .*")
file(STRINGS "${TAMPI_INCLUDE_DIR}/TAMPI_Decl.h" version-minor-line
  REGEX "#define TAMPI_VERSION_MINOR .*")

if (NOT version-major-line OR NOT version-minor-line)
  message(FATAL "Cannot find TAMPI version lines")
endif()

string(REGEX REPLACE "^#define[ \t]+TAMPI_VERSION_MAJOR[ \t]+([0-9]+)$" "\\1"
  TAMPI_VERSION_MAJOR ${version-major-line})
string(REGEX REPLACE "^#define[ \t]+TAMPI_VERSION_MINOR[ \t]+([0-9]+)$" "\\1"
  TAMPI_VERSION_MINOR ${version-minor-line})

set(TAMPI_VERSION ${TAMPI_VERSION_MAJOR}.${TAMPI_VERSION_MINOR} CACHE STRING "TAMPI Version")

message(STATUS "Found TAMPI version ${TAMPI_VERSION}")

add_library(Tampi::tampi-c SHARED IMPORTED)
set_target_properties(Tampi::tampi-c PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${TAMPI_INCLUDE_DIR}"
  IMPORTED_LOCATION ${TAMPI_LIBRARY})
