# Copyright (c) 2022-2024 Barcelona Supercomputing Center (BSC)
# SPDX-License-Identifier: GPL-3.0-or-later

set(OPENMPV_FLAG "-fopenmp=libompv")

include(CheckCCompilerFlag)

set(CMAKE_REQUIRED_LINK_OPTIONS "${OPENMPV_FLAG}")
check_c_compiler_flag("${OPENMPV_FLAG}" OPENMPV_FOUND)

if(NOT OPENMPV_FOUND)
  message(STATUS "Compiler doesn't support -fopenmp=libompv")
  return()
endif()

if(NOT TARGET OpenMPV)
  add_library(OpenMPV INTERFACE)
  target_compile_options(OpenMPV INTERFACE "${OPENMPV_FLAG}")
  target_link_options(OpenMPV INTERFACE "${OPENMPV_FLAG}")
endif()
