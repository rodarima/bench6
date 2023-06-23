# Copyright (c) 2021-2023 Barcelona Supercomputing Center (BSC)
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.20)

project(BENCH6 LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

add_compile_options(
  -Wall
  -Wextra
  -Wformat
  -Wmissing-prototypes
  -Wstrict-prototypes
  -Wold-style-definition
  -pedantic
  -Werror
)

set(CMAKE_C_VISIBILITY_PRESET hidden)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED TRUE)
set(CMAKE_C_EXTENSIONS FALSE)

add_definitions(-D_POSIX_C_SOURCE=200809L)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel Asan UBsan." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

include(GNUInstallDirs)

# Extra build type for AddressSanitizer (Asan)
set(CMAKE_C_FLAGS_ASAN "${CMAKE_C_FLAGS_DEBUG} \
  -fsanitize=address \
  -fno-optimize-sibling-calls \
  -fsanitize-address-use-after-scope \
  -fno-omit-frame-pointer"
  CACHE STRING "Flags used by the C compiler during AddressSanitizer builds." FORCE)

set(CMAKE_C_FLAGS_UBSAN "${CMAKE_C_FLAGS_DEBUG} -fsanitize=undefined"
  CACHE STRING "Flags used by the C compiler during UndefinedBehaviorSanitizer builds." FORCE)

# Required for clock_gettime() in glibc <= 2.17
include(CheckLibraryExists)
check_library_exists(c clock_gettime "" HAVE_CLOCK_GETTIME)

if (NOT HAVE_CLOCK_GETTIME)
  check_library_exists(rt clock_gettime "time.h" HAVE_CLOCK_GETTIME_RT)
  if (HAVE_CLOCK_GETTIME_RT)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lrt")
    set(CMAKE_C_STANDARD_LIBRARIES "${CMAKE_C_STANDARD_LIBRARIES} -lrt")
  else()
    message(FATAL_ERROR "cannot find clock_gettime()")
  endif()
endif()

set(USE_SIMD ON CACHE BOOL "Define SIMD and build with -fopenmp-simd")
if(USE_SIMD)
  add_definitions(-DSIMD)
  add_compile_options(-fopenmp-simd)
endif()

# Use <PackageName>_ROOT variables if available, commonly used by MPI
# installations
if(POLICY CMP0074)
cmake_policy(SET CMP0074 NEW)
endif()
find_package(MPI)
find_package(Nanos6)
find_package(Nodes)
find_package(Tampi)

set_property(GLOBAL PROPERTY bench6_list "")

macro(mk_bench NAME)
  if(NOT "${NAME}" MATCHES "b6_.*")
    message(FATAL_ERROR "benchmark name must begin with b6_: ${NAME}")
  endif()
  add_executable(${NAME})

  get_property(BENCH6_LIST GLOBAL PROPERTY bench6_list)
  message(STATUS "Before BENCH6_LIST=${BENCH6_LIST}")
  list(APPEND BENCH6_LIST ${NAME})
  message(STATUS "After BENCH6_LIST=${BENCH6_LIST}")
  set_property(GLOBAL PROPERTY bench6_list "${BENCH6_LIST}")

  install(TARGETS ${NAME} RUNTIME DESTINATION bin)
endmacro()

add_subdirectory(src)

include(FeatureSummary)
feature_summary(WHAT ALL)
