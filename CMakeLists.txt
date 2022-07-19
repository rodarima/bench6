# Copyright (c) 2022 Barcelona Supercomputing Center (BSC)
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.10)

project(BENCH6 LANGUAGES C)

add_compile_options(-Wall -Wextra -Wformat
  -Wmissing-prototypes -Wstrict-prototypes
  #-Wconversion -Wsign-conversion
  -Wold-style-definition -pedantic
  -Werror
)

set(CMAKE_C_COMPILER "clang")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_EXTENSIONS FALSE)

add_executable(bench6
  src/main.c
  src/common.c
  src/sched_get.c
  src/sched_add.c
  src/register_deps.c)

target_compile_options(bench6 PRIVATE -fompss-2)
target_link_options(bench6 PRIVATE -fompss-2)

install(TARGETS bench6 RUNTIME DESTINATION bin)
