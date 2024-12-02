# -*- Makefile -*-

# - Compile time options -----------------------------------------------
#
# -DHPCG_CONTIGUOUS_ARRAYS Define to have sparse matrix arrays long and contiguous
# -DHPCG_DEBUG       	   Define to enable debugging output
# -DHPCG_DETAILED_DEBUG    Define to enable very detailed debugging output
#
# By default HPCG will:
#    *) Not generate debugging output.
#
HPCG_OPTS     =

CXX          = clang++
CXXFLAGS     = $(HPCG_OPTS) -O3 -mavx -fopenmp -I$(TAMPI_HOME)/include
MPICXX       = mpiicpc
LINKER       = $(WRAPPERS) $(WRAPPERS) $(MPICXX)
LINKFLAGS    = -L$(TAMPI_HOME)/lib -l:libtampi.a

# MPI Wrappers
WRAPPERS=I_MPI_CXX=$(CXX) MPICH_CXX=$(CXX) OMPI_CXX=$(CXX)

HPCG_DEPS = src/CG.o \
	    src/CG_ref.o \
	    src/TestCG.o \
	    src/ComputeResidual.o \
	    src/ExchangeHalo_ref.o \
	    src/ExchangeHalo.o \
	    src/GenerateGeometry.o \
	    src/GenerateProblem.o \
	    src/GenerateProblem_ref.o \
	    src/CheckProblem.o \
	    src/MixedBaseCounter.o \
	    src/OptimizeProblem.o \
	    src/ReadHpcgDat.o \
	    src/ReportResults.o \
	    src/SetupHalo.o \
	    src/SetupHalo_ref.o \
	    src/TestSymmetry.o \
	    src/TestNorms.o \
	    src/WriteProblem.o \
	    src/YAML_Doc.o \
	    src/YAML_Element.o \
	    src/ComputeDotProduct.o \
	    src/ComputeDotProduct_ref.o \
	    src/mytimer.o \
	    src/ComputeOptimalShapeXYZ.o \
	    src/ComputeSPMV.o \
	    src/ComputeSPMV_ref.o \
	    src/ComputeSYMGS.o \
	    src/ComputeSYMGS_ref.o \
	    src/ComputeWAXPBY.o \
	    src/ComputeWAXPBY_ref.o \
	    src/ComputeMG_ref.o \
	    src/ComputeMG.o \
	    src/ComputeProlongation_ref.o \
	    src/ComputeRestriction_ref.o \
	    src/ComputeProlongation.o \
	    src/ComputeRestriction.o \
	    src/CheckAspectRatio.o \
	    src/OutputFile.o \
	    src/GenerateCoarseProblem.o \
	    src/init.o \
	    src/finalize.o \
	    src/Serialization.o

# These header files are included in many source files, so we recompile every file if one or more of these header is modified.
PRIMARY_HEADERS = src/Geometry.hpp src/SparseMatrix.hpp src/Vector.hpp src/CGData.hpp \
                  src/MGData.hpp src/hpcg.hpp

all: bin/xhpcg

bin/xhpcg: src/main.o $(HPCG_DEPS)
	$(LINKER) $(CXXFLAGS) src/main.o $(HPCG_DEPS) $(HPCG_LIBS) $(LINKFLAGS) -o bin/xhpcg

clean:
	rm -f src/*.o bin/xhpcg

.PHONY: all clean

src/main.o: src/main.cpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/CG.o: src/CG.cpp src/CG.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/CG_ref.o: src/CG_ref.cpp src/CG_ref.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/TestCG.o: src/TestCG.cpp src/TestCG.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeResidual.o: src/ComputeResidual.cpp src/ComputeResidual.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ExchangeHalo_ref.o: src/ExchangeHalo_ref.cpp src/ExchangeHalo_ref.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ExchangeHalo.o: src/ExchangeHalo.cpp src/ExchangeHalo.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/GenerateGeometry.o: src/GenerateGeometry.cpp src/GenerateGeometry.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/GenerateProblem.o: src/GenerateProblem.cpp src/GenerateProblem.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/GenerateProblem_ref.o: src/GenerateProblem_ref.cpp src/GenerateProblem_ref.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/CheckProblem.o: src/CheckProblem.cpp src/CheckProblem.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/MixedBaseCounter.o: src/MixedBaseCounter.cpp src/MixedBaseCounter.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/OptimizeProblem.o: src/OptimizeProblem.cpp src/OptimizeProblem.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ReadHpcgDat.o: src/ReadHpcgDat.cpp src/ReadHpcgDat.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ReportResults.o: src/ReportResults.cpp src/ReportResults.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/SetupHalo.o: src/SetupHalo.cpp src/SetupHalo.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/SetupHalo_ref.o: src/SetupHalo_ref.cpp src/SetupHalo_ref.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/TestSymmetry.o: src/TestSymmetry.cpp src/TestSymmetry.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/TestNorms.o: src/TestNorms.cpp src/TestNorms.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/WriteProblem.o: src/WriteProblem.cpp src/WriteProblem.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/YAML_Doc.o: src/YAML_Doc.cpp src/YAML_Doc.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/YAML_Element.o: src/YAML_Element.cpp src/YAML_Element.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeDotProduct.o: src/ComputeDotProduct.cpp src/ComputeDotProduct.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeDotProduct_ref.o: src/ComputeDotProduct_ref.cpp src/ComputeDotProduct_ref.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/finalize.o: src/finalize.cpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/init.o: src/init.cpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/mytimer.o: src/mytimer.cpp src/mytimer.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeOptimalShapeXYZ.o: src/ComputeOptimalShapeXYZ.cpp src/ComputeOptimalShapeXYZ.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeSPMV.o: src/ComputeSPMV.cpp src/ComputeSPMV.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeSPMV_ref.o: src/ComputeSPMV_ref.cpp src/ComputeSPMV_ref.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeSYMGS.o: src/ComputeSYMGS.cpp src/ComputeSYMGS.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeSYMGS_ref.o: src/ComputeSYMGS_ref.cpp src/ComputeSYMGS_ref.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeWAXPBY.o: src/ComputeWAXPBY.cpp src/ComputeWAXPBY.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeWAXPBY_ref.o: src/ComputeWAXPBY_ref.cpp src/ComputeWAXPBY_ref.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeMG_ref.o: src/ComputeMG_ref.cpp src/ComputeMG_ref.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeMG.o: src/ComputeMG.cpp src/ComputeMG.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeProlongation_ref.o: src/ComputeProlongation_ref.cpp src/ComputeProlongation_ref.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeRestriction_ref.o: src/ComputeRestriction_ref.cpp src/ComputeRestriction_ref.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeProlongation.o: src/ComputeProlongation.cpp src/ComputeProlongation.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/ComputeRestriction.o: src/ComputeRestriction.cpp src/ComputeRestriction.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/GenerateCoarseProblem.o: src/GenerateCoarseProblem.cpp src/GenerateCoarseProblem.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/CheckAspectRatio.o: src/CheckAspectRatio.cpp src/CheckAspectRatio.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/OutputFile.o: src/OutputFile.cpp src/OutputFile.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

src/Serialization.o: src/Serialization.cpp src/Serialization.hpp $(PRIMARY_HEADERS)
	$(WRAPPERS) $(MPICXX) -c $(CXXFLAGS) -Isrc $< -o $@

