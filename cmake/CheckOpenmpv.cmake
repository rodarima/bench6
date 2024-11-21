# Copyright (c) 2022-2023 Barcelona Supercomputing Center (BSC)
# SPDX-License-Identifier: GPL-3.0-or-later

include(CheckCCompilerFlag)

set(CMAKE_REQUIRED_LINK_OPTIONS "-fopenmp=libompv")
check_c_compiler_flag("-fopenmp=libompv" OPENMPV_COMPILER_FOUND)
