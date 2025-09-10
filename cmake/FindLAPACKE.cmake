include(GNUInstallDirs)

if (TARGET MKL::MKL)
  # If MKL is found, create an alias for CBLAS to MKL::MKL
  target_compile_definitions(MKL::MKL INTERFACE USE_MKL)
  add_library(LAPACKE ALIAS MKL::MKL)
  set(LAPACKE_FOUND TRUE)
  return()
endif()

find_library(LAPACKE_LIBRARY NAMES lapacke)
find_path(LAPACKE_INCLUDE_DIR lapacke.h)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(LAPACKE DEFAULT_MSG
  LAPACKE_LIBRARY LAPACKE_INCLUDE_DIR)

if(NOT LAPACKE_FOUND)
  return()
endif()

if(TARGET LAPACKE)
  return()
endif()

add_library(LAPACKE SHARED IMPORTED)
set_target_properties(LAPACKE PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${LAPACKE_INCLUDE_DIR}"
  IMPORTED_LOCATION ${LAPACKE_LIBRARY})
