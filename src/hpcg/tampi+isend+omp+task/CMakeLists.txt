
set(SOURCES src/main.cpp src/CG.cpp src/CG_ref.cpp src/TestCG.cpp
  src/ComputeResidual.cpp src/ExchangeHalo.cpp src/ExchangeHalo_ref.cpp src/GenerateGeometry.cpp
  src/GenerateProblem.cpp src/GenerateProblem_ref.cpp src/CheckProblem.cpp
  src/OptimizeProblem.cpp src/ReadHpcgDat.cpp src/ReportResults.cpp
  src/SetupHalo.cpp src/SetupHalo_ref.cpp src/TestSymmetry.cpp src/TestNorms.cpp
  src/WriteProblem.cpp src/YAML_Doc.cpp src/YAML_Element.cpp
  src/ComputeDotProduct.cpp src/ComputeDotProduct_ref.cpp src/finalize.cpp
  src/init.cpp src/mytimer.cpp src/ComputeSPMV.cpp src/ComputeSPMV_ref.cpp
  src/ComputeSYMGS.cpp src/ComputeSYMGS_ref.cpp src/ComputeWAXPBY.cpp
  src/ComputeWAXPBY_ref.cpp src/ComputeMG_ref.cpp src/ComputeMG.cpp
  src/ComputeProlongation.cpp src/ComputeProlongation_ref.cpp
  src/ComputeRestriction.cpp src/ComputeRestriction_ref.cpp
  src/GenerateCoarseProblem.cpp src/ComputeOptimalShapeXYZ.cpp
  src/MixedBaseCounter.cpp src/CheckAspectRatio.cpp src/OutputFile.cpp
  src/Serialization.cpp)

# CXXFLAGS     = $(HPCG_OPTS) -O3 -mavx -fopenmp -I$(TAMPI_HOME)/include

if(MPI_FOUND AND TAMPI_FOUND)
  if(OPENMPV_FOUND)
    mk_bench(b6_hpcg_tampi_ompv)
    target_sources(b6_hpcg_tampi_ompv PRIVATE ${SOURCES})
    target_link_libraries(b6_hpcg_tampi_ompv
      PRIVATE MPI::MPI_C Tampi::tampi-c OpenMPV)
  endif()
endif()
