
# @HEADER
# ***********************************************************************
#  
#                HPCCG: Simple Conjugate Gradient Benchmark Code
#                  Copyright (2006) Sandia Corporation
#  
#  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
#  license for use of this work by or on behalf of the U.S. Government.
#  
#  BSD 3-Clause License
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#  
#  * Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#  
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  
#  * Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#  
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  
#  Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
#  
#  ************************************************************************
# @HEADER

# Simple hand-tuned makefile.  Modify as necessary for your environment.
# Questions? Contact Mike Heroux (maherou@sandia.gov).
#

#
# 0) Specify compiler and linker:

CXX=icpc

#IA32 with GCC: 
ICPC_VEC_FLAGS=-Wall
CPP_OPT_FLAGS = -O3 -march=native -mcpu=native -mtune=native $(ICPC_VEC_FLAGS)

#
# 1) OpenMP Compiler argument
#    GCC and Intel compilers require -fopenmp and -openmp, resp.  Other compilers may differ.

OMP_FLAGS = -fopenmp

# 7) System libraries: (May need to add -lg2c before -lm)

SYS_LIB =-lm

#
# 6) Specify name if executable (optional):

TARGET = test_HPCCG

################### Derived Quantities (no modification required) ##############

CPPFLAGS = -Iinclude -DWALL
CXXFLAGS= $(CPP_OPT_FLAGS)

LIB_PATHS= $(SYS_LIB)

INC_PATHS= include

OBJ_DIR := obj

SRC_OMP_DIR := src/mpi_interop_ompss
SRC_OMP_FILES := $(wildcard $(SRC_OMP_DIR)/*.cpp)
OMP_FILES := $(patsubst $(SRC_OMP_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_OMP_FILES))
COMMON_DIR := src/common
COMMON_FILES := $(wildcard $(COMMON_DIR)/*.cpp)
COMMON_OMP_FILES := $(patsubst $(COMMON_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(COMMON_FILES))

.PHONY: all
all: $(TARGET)

# OMP OBJS
$(OBJ_DIR)/%.o: $(SRC_OMP_DIR)/%.cpp
	   $(CXX) $(CPPFLAGS) $(USE_MPI) $(OMP_FLAGS) $(CXXFLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: $(COMMON_DIR)/%.cpp
	   $(CXX) $(CPPFLAGS) $(USE_MPI) $(OMP_FLAGS) $(CXXFLAGS) -c -o $@ $<

# TARGETS 
$(TARGET): $(COMMON_OMP_FILES) $(OMP_FILES)
	$(CXX) $(COMMON_OMP_FILES) $(OMP_FILES) $(OMP_FLAGS) $(LIB_PATHS) -o $@
test:
	@echo "Not implemented yet..."

.PHONY: clean
clean:
	@rm -f obj/* *~ $(TARGET)
