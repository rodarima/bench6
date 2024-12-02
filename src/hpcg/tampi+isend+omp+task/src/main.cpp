
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/*!
 @file main.cpp

 HPCG routine
 */

// Main routine of a program that calls the HPCG conjugate gradient
// solver to solve the problem, and then prints results.

#include <mpi.h>

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#ifdef HPCG_DETAILED_DEBUG
using std::cin;
#endif
using std::endl;
using std::isnan;

#include <vector>
#include <sstream>

#include "hpcg.hpp"

#include "CheckAspectRatio.hpp"
#include "GenerateGeometry.hpp"
#include "GenerateProblem.hpp"
#include "GenerateCoarseProblem.hpp"
#include "SetupHalo.hpp"
#include "CheckProblem.hpp"
#include "ExchangeHalo.hpp"
#include "OptimizeProblem.hpp"
#include "WriteProblem.hpp"
#include "ReportResults.hpp"
#include "mytimer.hpp"
#include "ComputeSPMV_ref.hpp"
#include "ComputeMG_ref.hpp"
#include "ComputeResidual.hpp"
#include "CG.hpp"
#include "CG_ref.hpp"
#include "Geometry.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"
#include "CGData.hpp"
#include "TestCG.hpp"
#include "TestSymmetry.hpp"
#include "TestNorms.hpp"
#include "Serialization.hpp"

/*!
  Main driver program: Construct synthetic problem, run V&V tests, compute benchmark parameters, run benchmark, report results.

  @param[in]  argc Standard argument count.  Should equal 1 (no arguments passed in) or 4 (nx, ny, nz passed in)
  @param[in]  argv Standard argument array.  If argc==1, argv is unused.  If argc==4, argv[1], argv[2], argv[3] will be interpreted as nx, ny, nz, resp.

  @return Returns zero on success and a non-zero value otherwise.

*/

const int MaxIters = 50;

// Use this array for collecting timing information
std::vector< double > times(10,0.0);
const int numberOfMgLevels = 4; // Number of levels including first

static std::string createPathStr(const char *dir, const char *fname, int rank, int level) {
  std::stringstream ss;
  ss << dir << "/" << fname << "." << rank << "." << level << ".bin";
  return ss.str();
}

static int Init(
    HPCG_Params &params,
    Geometry *&geom, SparseMatrix &A, Vector &b, Vector &x, Vector &xexact,
    CGData &data, double &refTolerance) {

  // Measure the time to perform the whole init phase without serialization
  double global_init_time = mytimer();

  // Check if QuickPath option is enabled.
  // If the running time is set to zero, we minimize all paths through the program
  bool quickPath = (params.runningTime==0);

  int size = params.comm_size;
  int rank = params.comm_rank; // Number of MPI processes, My process ID

#ifdef HPCG_DETAILED_DEBUG
  if (size < 100 && rank==0) HPCG_fout << "Process "<<rank<<" of "<<size<<" is alive with " << params.numThreads << " threads." <<endl;

  if (rank==0) {
    char c;
    std::cout << "Press key to continue"<< std::endl;
    std::cin.get(c);
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  local_int_t nx,ny,nz;
  nx = (local_int_t)params.nx;
  ny = (local_int_t)params.ny;
  nz = (local_int_t)params.nz;
  int ierr = 0;  // Used to check return codes on function calls

  if (!params.noAspectRatio)
    ierr = CheckAspectRatio(0.125, nx, ny, nz, "local problem", rank==0);
  if (ierr)
    return ierr;

  /////////////////////////
  // Problem setup Phase //
  /////////////////////////

#ifdef HPCG_DEBUG
  double t1 = mytimer();
#endif

  // Construct the geometry and linear system
  geom = new Geometry;
  GenerateGeometry(size, rank, params.numThreads, params.pz, params.zl, params.zu, nx, ny, nz, params.npx, params.npy, params.npz, geom);

  if (!params.noAspectRatio)
    ierr = CheckAspectRatio(0.125, geom->npx, geom->npy, geom->npz, "process grid", rank==0);
  if (ierr)
    return ierr;

  double setup_time = mytimer();

  InitializeSparseMatrix(A, geom);

  GenerateProblem(A, &b, &x, &xexact);
  SetupHalo(A);
  SparseMatrix * curLevelMatrix = &A;
  for (int level = 1; level< numberOfMgLevels; ++level) {
    GenerateCoarseProblem(*curLevelMatrix);
    curLevelMatrix = curLevelMatrix->Ac; // Make the just-constructed coarse grid the next level
  }

  setup_time = mytimer() - setup_time; // Capture total time of setup
  times[9] = setup_time; // Save it for reporting

  curLevelMatrix = &A;
  Vector * curb = &b;
  Vector * curx = &x;
  Vector * curxexact = &xexact;
  for (int level = 0; level< numberOfMgLevels; ++level) {
     CheckProblem(*curLevelMatrix, curb, curx, curxexact);
     curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
     curb = 0; // No vectors after the top level
     curx = 0;
     curxexact = 0;
  }

  InitializeSparseCGData(A, data);

  ////////////////////////////////////
  // Reference SpMV+MG Timing Phase //
  ////////////////////////////////////

  // Call Reference SpMV and MG. Compute Optimization time as ratio of times in these routines

  local_int_t nrow = A.localNumberOfRows;
  local_int_t ncol = A.localNumberOfColumns;

  Vector x_overlap, b_computed;
  InitializeVector(x_overlap, ncol); // Overlapped copy of x vector
  InitializeVector(b_computed, nrow); // Computed RHS vector

  // Record execution time of reference SpMV and MG kernels for reporting times
  // First load vector with random values
  FillRandomVector(x_overlap);

  int numberOfCalls = 10;
  if (quickPath) numberOfCalls = 1; //QuickPath means we do on one call of each block of repetitive code
  double t_begin = mytimer();
  for (int i=0; i< numberOfCalls; ++i) {
    ierr = ComputeSPMV_ref(A, x_overlap, b_computed); // b_computed = A*x_overlap
    if (ierr) HPCG_fout << "Error in call to SpMV: " << ierr << ".\n" << endl;
    ierr = ComputeMG_ref(A, b_computed, x_overlap); // b_computed = Minv*y_overlap
    if (ierr) HPCG_fout << "Error in call to MG: " << ierr << ".\n" << endl;
  }
  times[8] = (mytimer() - t_begin)/((double) numberOfCalls);  // Total time divided by number of calls.
#ifdef HPCG_DEBUG
  if (rank==0) HPCG_fout << "Total SpMV+MG timing phase execution time in main (sec) = " << mytimer() - t1 << endl;
#endif

  ///////////////////////////////
  // Reference CG Timing Phase //
  ///////////////////////////////

#ifdef HPCG_DEBUG
  t1 = mytimer();
#endif

  int niters = 0;
  double normr = 0.0;
  double normr0 = 0.0;
  int refMaxIters = MaxIters;
  numberOfCalls = 1; // Only need to run the residual reduction analysis once

  // Compute the residual reduction for the natural ordering and reference kernels
  std::vector< double > ref_times(9,0.0);
  double tolerance = 0.0; // Set tolerance to zero to make all runs do maxIters iterations
  int err_count = 0;
  for (int i=0; i< numberOfCalls; ++i) {
    ZeroVector(x);
    ierr = CG_ref( A, data, b, x, refMaxIters, tolerance, niters, normr, normr0, &ref_times[0], true);
    if (ierr) ++err_count; // count the number of errors in CG
  }
  if (rank == 0 && err_count) HPCG_fout << err_count << " error(s) in call(s) to reference CG." << endl;
  refTolerance = normr / normr0;

  global_init_time = mytimer() - global_init_time;
  if (params.store_path) {
    // Print only when running init phase only
    if (rank == 0)
      printf("time %e\n", global_init_time);

    ParamsSerialize(params, createPathStr(params.store_path, "HPCG_Params", rank, 0));
    ToleranceSerialize(refTolerance, createPathStr(params.store_path, "refTolerance", rank, 0));

    // Clean up
    DeleteMatrix(A); // This delete will recursively delete all coarse grid data
    DeleteVector(x);
    DeleteVector(b);
    DeleteVector(xexact);
  }
  DeleteVector(x_overlap);
  DeleteVector(b_computed);

  return 0;

}

static int Run(
    HPCG_Params &params,
    Geometry *&geom, SparseMatrix &A, Vector &b, Vector &x, Vector &xexact,
    CGData &data, double &refTolerance) {

  int size = params.comm_size;
  int rank = params.comm_rank; // Number of MPI processes, My process ID

  double deserial_time = 0.0;
  if (params.load_path) {
    deserial_time = mytimer();
    HPCG_Params deserial_params;
    ParamsDeserialize(deserial_params, createPathStr(params.load_path, "HPCG_Params", rank, 0));
    ToleranceDeserialize(refTolerance, createPathStr(params.load_path, "refTolerance", rank, 0));
    deserial_time = mytimer() - deserial_time;
    if (rank == 0)
      printf("deserial_time %e\n", deserial_time);

    // Check we are not deserializing other problem size
    die_if(deserial_params.comm_size != params.comm_size, "deserial params.comm_size does not match");
    die_if(deserial_params.comm_rank != params.comm_rank, "deserial params.comm_rank does not match");
    // die_if(deserial_params.numThreads != params.numThreads, "deserial params.numThreads does not match");
    die_if(deserial_params.nx != params.nx, "deserial params.nx does not match");
    die_if(deserial_params.ny != params.ny, "deserial params.ny does not match");
    die_if(deserial_params.nz != params.nz, "deserial params.nz does not match");
    // die_if(deserial_params.runningTime != params.runningTime, "deserial params.runningTime does not match");
    die_if(deserial_params.npx != params.npx, "deserial params.npx does not match");
    die_if(deserial_params.npy != params.npy, "deserial params.npy does not match");
    die_if(deserial_params.npz != params.npz, "deserial params.npz does not match");
    die_if(deserial_params.pz != params.pz, "deserial params.pz does not match");
    die_if(deserial_params.zl != params.zl, "deserial params.zl does not match");
    die_if(deserial_params.zu != params.zu, "deserial params.zu does not match");
    // Do not check numBlocks since it does not have huge
    // impact in the result apart from floating reductions
    // die_if(deserial_params.numBlocks != params.numBlocks, "deserial params.numBlocks does not match");
    // die_if(deserial_params.numNeighComms != params.numNeighComms, "deserial params.numNeighComms does not match");
    // die_if(deserial_params.noAspectRatio != params.noAspectRatio, "deserial params.noAspectRatio does not match");
    DeleteParamsData(deserial_params);
  }
  if (params.load_path || !isnan(params.tolerance)) {
    double init_time = mytimer();
    local_int_t nx,ny,nz;
    nx = (local_int_t)params.nx;
    ny = (local_int_t)params.ny;
    nz = (local_int_t)params.nz;

    geom = new Geometry;
    GenerateGeometry(size, rank, params.numThreads, params.pz, params.zl, params.zu, nx, ny, nz, params.npx, params.npy, params.npz, geom);

    InitializeSparseMatrix(A, geom);

    GenerateProblem(A, &b, &x, &xexact);
    SetupHalo(A);
    SparseMatrix *curLevelMatrix = &A;
    for (int level = 1; level < numberOfMgLevels; ++level) {
      GenerateCoarseProblem(*curLevelMatrix);
      curLevelMatrix = curLevelMatrix->Ac; // Make the just-constructed coarse grid the next level
    }

    InitializeSparseCGData(A, data);
    init_time = mytimer() - init_time;
    if (rank == 0) {
      printf("init_time %e\n", init_time);
      printf("init_time+deserial_time %e\n", init_time+deserial_time);
    }
  }

#ifdef HPCG_DEBUG
  double t1 = mytimer();
#endif

  bool quickPath = (params.runningTime==0);
  int ierr = 0;  // Used to check return codes on function calls

  // Call user-tunable set up function.
  double t7 = mytimer();
  OptimizeProblem(A, data, b, x, xexact, params.numBlocks, params.numNeighComms);
  t7 = mytimer() - t7;
  times[7] = t7;
#ifdef HPCG_DEBUG
  if (rank==0) HPCG_fout << "Total problem setup time in main (sec) = " << mytimer() - t1 << endl;
#endif

#ifdef HPCG_DETAILED_DEBUG
  if (geom->size == 1) WriteProblem(*geom, A, b, x, xexact);
#endif

  //////////////////////////////
  // Validation Testing Phase //
  //////////////////////////////

#ifdef HPCG_DEBUG
  t1 = mytimer();
#endif
  TestCGData testcg_data;
  testcg_data.count_pass = testcg_data.count_fail = 0;
  TestCG(A, data, b, x, testcg_data, params.numNeighComms);

  TestSymmetryData testsymmetry_data;
  TestSymmetry(A, b, xexact, testsymmetry_data, params.numNeighComms);

#ifdef HPCG_DEBUG
  if (rank==0) HPCG_fout << "Total validation (TestCG and TestSymmetry) execution time in main (sec) = " << mytimer() - t1 << endl;
#endif

#ifdef HPCG_DEBUG
  t1 = mytimer();
#endif

  //////////////////////////////
  // Optimized CG Setup Phase //
  //////////////////////////////

  int niters = 0;
  double normr = 0.0;
  double normr0 = 0.0;
  int err_count = 0;
  int tolerance_failures = 0;
  int global_failure = 0; // assume all is well: no failures

  int refMaxIters = MaxIters;
  int optMaxIters = 10*refMaxIters;
  int optNiters = refMaxIters;
  double opt_worst_time = 0.0;

  int numberOfCalls = 1; // Only need to run the residual reduction analysis once

  std::vector< double > opt_times(9,0.0);

  // Compute the residual reduction and residual count for the user ordering and optimized kernels.
  for (int i=0; i< numberOfCalls; ++i) {
    ZeroVector(x); // start x at all zeros
    double last_cummulative_time = opt_times[0];
    ierr = CG( A, data, b, x, optMaxIters, refTolerance, niters, normr, normr0, &opt_times[0], params.numNeighComms, true);
    if (ierr) ++err_count; // count the number of errors in CG
#ifdef HPCG_DEBUG
    if (rank == 0)
      std::cout << "normr " << normr << " normr0 " << normr0 << " refTolerance " << refTolerance << "\n";
#endif
    if (normr / normr0 > refTolerance) ++tolerance_failures; // the number of failures to reduce residual

    // pick the largest number of iterations to guarantee convergence
    if (niters > optNiters) optNiters = niters;

    double current_time = opt_times[0] - last_cummulative_time;
    if (current_time > opt_worst_time) opt_worst_time = current_time;
  }

// Get the absolute worst time across all MPI ranks (time in CG can be different)
  double local_opt_worst_time = opt_worst_time;
  MPI_Allreduce(&local_opt_worst_time, &opt_worst_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if (rank == 0 && err_count) HPCG_fout << err_count << " error(s) in call(s) to optimized CG." << endl;
  if (tolerance_failures) {
    global_failure = 1;
    if (rank == 0)
      HPCG_fout << "Failed to reduce the residual " << tolerance_failures << " times." << endl;
  }

  ///////////////////////////////
  // Optimized CG Timing Phase //
  ///////////////////////////////

  // Here we finally run the benchmark phase
  // The variable total_runtime is the target benchmark execution time in seconds

  double total_runtime = params.runningTime;
  int numberOfCgSets = int(total_runtime / opt_worst_time) + 1; // Run at least once, account for rounding

#ifdef HPCG_DEBUG
  if (rank==0) {
    HPCG_fout << "Projected running time: " << total_runtime << " seconds" << endl;
    HPCG_fout << "Number of CG sets: " << numberOfCgSets << endl;
  }
#endif

  /* This is the timed run for a specified amount of time. */

  optMaxIters = optNiters;
  double optTolerance = 0.0;  // Force optMaxIters iterations
  TestNormsData testnorms_data;
  testnorms_data.samples = numberOfCgSets;
  testnorms_data.values = new double[numberOfCgSets];

  for (int i=0; i< numberOfCgSets; ++i) {
    ZeroVector(x); // Zero out x
    ierr = CG( A, data, b, x, optMaxIters, optTolerance, niters, normr, normr0, &times[0], params.numNeighComms, true);
    if (ierr) HPCG_fout << "Error in call to CG: " << ierr << ".\n" << endl;
    if (rank==0) HPCG_fout << "Call [" << i << "] Scaled Residual [" << normr/normr0 << "]" << endl;
    testnorms_data.values[i] = normr/normr0; // Record scaled residual from this run
  }

  // Compute difference between known exact solution and computed solution
  // All processors are needed here.
#ifdef HPCG_DEBUG
  double residual = 0;
  ierr = ComputeResidual(A.localNumberOfRows, x, xexact, residual);
  if (ierr) HPCG_fout << "Error in call to compute_residual: " << ierr << ".\n" << endl;
  if (rank==0) HPCG_fout << "Difference between computed and exact  = " << residual << ".\n" << endl;
#endif

  // Test Norm Results
  ierr = TestNorms(testnorms_data);

  ////////////////////
  // Report Results //
  ////////////////////

  // Report results to YAML file
  ReportResults(A, numberOfMgLevels, numberOfCgSets, refMaxIters, optMaxIters, &times[0], testcg_data, testsymmetry_data, testnorms_data, global_failure, quickPath);

  // Clean up
  DeleteMatrix(A); // This delete will recursively delete all coarse grid data
  DeleteVector(x);
  DeleteVector(b);
  DeleteVector(xexact);
  delete [] testnorms_data.values;

  HPCG_Finalize();

  return 0;
}

int main(int argc, char * argv[]) {
  HPCG_Params params;
  Geometry *geom = 0;
  SparseMatrix A;
  Vector b, x, xexact;
  CGData data;
  double refTolerance = 0.0;

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  assert(provided == MPI_THREAD_MULTIPLE);

  HPCG_Init(&argc, &argv, params);

  die_if((params.load_path && !isnan(params.tolerance)), "Cannot use both at the same time");
  die_if((params.store_path && !isnan(params.tolerance)), "Cannot use both at the same time");
  die_if((params.load_path && params.store_path), "Cannot use both at the same time");


  if (!params.load_path && isnan(params.tolerance)) {
    int ierr = Init(params, geom, A, b, x, xexact, data, refTolerance);
    die_if(ierr, "Something were wrong in initialization phase");
  }

  if (!params.store_path) {
    double tmpTolerance = isnan(params.tolerance) ? refTolerance : params.tolerance;
    int ierr = Run(params, geom, A, b, x, xexact, data, tmpTolerance);
    die_if(ierr, "Something were wrong in run phase");
  }

  DeleteCGData(data);
  DeleteParamsData(params);

  // Finish up
  MPI_Finalize();
  return 0;
}
