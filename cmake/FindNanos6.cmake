# Copyright (c) 2024 Barcelona Supercomputing Center (BSC)
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Searchs for libnanos6 and checks if the -fompss-2=libnanos6 flag is supported
# by the compiler.
#
# Sets the variable NANOS6_FOUND when BOTH libnanos6 is found and the compiler
# supports the -fompss-2=libnanos6 flag.
#
# The target Nanos6::nanos6 is defined when both checks are passed, and includes
# a rule to add the compile and link time flags.

include(GNUInstallDirs)

set(NANOS6_FLAG "-fompss-2=libnanos6")

if(DEFINED ENV{NANOS6_HOME})
  set(NANOS6_HOME "$ENV{NANOS6_HOME}")

  # Ensure the compiler supports libnanos6
  include(CheckCCompilerFlag)

  # Also set the linker flags, as otherwise the check will fail
  # due to undefined symbols in the final program.
  set(CMAKE_REQUIRED_LINK_OPTIONS "${NANOS6_FLAG}")
  check_c_compiler_flag("${NANOS6_FLAG}" NANOS6_FLAG_SUPPORTED)

  find_library(NANOS6_LIBRARY NAMES nanos6 PATHS "${NANOS6_HOME}/lib" NO_DEFAULT_PATH)
  find_path(NANOS6_INCLUDE_DIR nanos6.h PATHS "${NANOS6_HOME}/include" NO_DEFAULT_PATH)
else()
  message(STATUS "NANOS6_HOME not set, refusing to search")
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Nanos6 DEFAULT_MSG
  NANOS6_HOME NANOS6_LIBRARY NANOS6_INCLUDE_DIR NANOS6_FLAG_SUPPORTED)

if(NOT NANOS6_FOUND)
  message(STATUS "Cannot find Nanos6 library")
  return()
endif()

if(NOT TARGET Nanos6::nanos6)
  add_library(Nanos6::nanos6 SHARED IMPORTED)
  set_target_properties(Nanos6::nanos6 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${NANOS6_INCLUDE_DIR}"
    IMPORTED_LOCATION ${NANOS6_LIBRARY})
  target_compile_options(Nanos6::nanos6 INTERFACE "${NANOS6_FLAG}")
  target_link_options(Nanos6::nanos6 INTERFACE "${NANOS6_FLAG}")
endif()
