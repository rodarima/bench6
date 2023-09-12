# Copyright (c) 2022-2023 Barcelona Supercomputing Center (BSC)
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Searchs for libnodes and checks if the -fompss-2=libnodes flag is supported in
# the compiler.
#
# Sets the variable NODES_FOUND when libnodes is found and NODES_FLAG_SUPPORTED
# when the -fompss-2=libnodes flag is supported.
#
# The target Nodes::nodes is defined when both checks are passed, and includes a
# rule to add the compile and link time flags.

include(GNUInstallDirs)

if(DEFINED ENV{NODES_HOME})
  set(NODES_HOME "$ENV{NODES_HOME}")
else()
  message(STATUS "NODES_HOME not set, refusing to search")
endif()

find_library(NODES_LIBRARY NAMES nodes PATHS "${NODES_HOME}/lib" NO_DEFAULT_PATH)
find_path(NODES_INCLUDE_DIR nodes.h PATHS "${NODES_HOME}/include" NO_DEFAULT_PATH)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Nodes DEFAULT_MSG
  NODES_LIBRARY NODES_INCLUDE_DIR)

if(NOT NODES_FOUND)
  message(STATUS "Cannot find NODES library")
  return()
endif()

# Also ensure the compiler supports libnodes
include(CheckCCompilerFlag)

set(NODES_FLAG "-fompss-2=libnodes")
# Also set the linker flags, as otherwise the check will fail due to undefined
# symbols in the final program.
set(CMAKE_REQUIRED_LINK_OPTIONS "${NODES_FLAG}")
check_c_compiler_flag("${NODES_FLAG}" NODES_FLAG_SUPPORTED)

if(NOT NODES_FLAG_SUPPORTED)
  message(STATUS "Compiler doesn't support ${NODES_FLAG} flag")
  return()
endif()

if(NOT TARGET Nodes::nodes)
  add_library(Nodes::nodes SHARED IMPORTED)
  set_target_properties(Nodes::nodes PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${NODES_INCLUDE_DIR}"
    IMPORTED_LOCATION ${NODES_LIBRARY})
  target_compile_options(Nodes::nodes INTERFACE "${NODES_FLAG}")
  target_link_options(Nodes::nodes INTERFACE "${NODES_FLAG}")
endif()
