
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
 @file TestSymmetry.cpp

 HPCG routine
 */

// The MPI include must be first for Windows platforms
#include <mpi.h>
#include <fstream>
#include <iostream>
#include <cfloat>
using std::endl;
#include <vector>
#include <cmath>

#include "hpcg.hpp"

#include "ComputeSPMV.hpp"
#include "ComputeMG.hpp"
#include "ComputeDotProduct.hpp"
#include "ComputeResidual.hpp"
#include "Geometry.hpp"
#include "SparseMatrix.hpp"
#include "TestSymmetry.hpp"

#include "ExchangeHalo.hpp"

/*!
  Tests symmetry-preserving properties of the sparse matrix vector multiply and multi-grid routines.

  @param[in]    geom   The description of the problem's geometry.
  @param[in]    A      The known system matrix
  @param[in]    b      The known right hand side vector
  @param[in]    xexact The exact solution vector
  @param[inout] testsymmetry_data The data structure with the results of the CG symmetry test including pass/fail information

  @return returns 0 upon success and non-zero otherwise

  @see ComputeDotProduct
  @see ComputeDotProduct_ref
  @see ComputeSPMV
  @see ComputeSPMV_ref
  @see ComputeMG
  @see ComputeMG_ref
*/
int TestSymmetry(SparseMatrix & A, Vector & b, Vector & xexact, TestSymmetryData & testsymmetry_data, int numNeighComms) {

 local_int_t ncol = A.localNumberOfColumns;

 Vector x_ncol, y_ncol, z_ncol;
 InitializeVector(x_ncol, ncol);
 InitializeVector(y_ncol, ncol);
 InitializeVector(z_ncol, ncol);

 // double t4 = 0.0; // Needed for dot-product call, otherwise unused
 testsymmetry_data.count_fail = 0;

 const BlockLayout & blkLayout = A.blkLayout;
 int numBlocks = blkLayout.numBlocks;
 const std::vector<BlockInfo> &blkInfo = blkLayout.blkInfo;

 // Test symmetry of matrix

 // First load vectors with random values
 FillRandomVector(x_ncol);
 FillRandomVector(y_ncol);

 double xNorm2, yNorm2;
 double xNorm2_local, yNorm2_local;
 double ANorm = 2 * 26.0;

 // Next, compute x'*A*y
 int ierr;
 for (int i = 0; i < numBlocks; i++) {
   size_t start = blkInfo[i].start;
   size_t end = blkInfo[i].end;
   ComputeDotProduct(start, end, y_ncol, y_ncol, yNorm2_local, A.isDotProductOptimized);
 }
 MPI_Allreduce(&yNorm2_local, &yNorm2, 1, MPI_DOUBLE, MPI_SUM,
     MPI_COMM_WORLD);
 yNorm2_local = 0.0;
 #pragma omp parallel
 #pragma omp single
  ExchangeHalo(A,y_ncol,numNeighComms);
 #pragma omp taskwait
 for (int i = 0; i < numBlocks; i++) {
   size_t start = blkInfo[i].start;
   size_t end = blkInfo[i].end;
   ierr = ComputeSPMV(A, start, end, y_ncol, z_ncol); // z_nrow = A*y_overlap
   if (ierr) HPCG_fout << "Error in call to SpMV: " << ierr << ".\n" << endl;
 }
 double xtAy = 0.0;
 double xtAy_local = 0.0;
 for (int i = 0; i < numBlocks; i++) {
   size_t start = blkInfo[i].start;
   size_t end = blkInfo[i].end;
   ierr = ComputeDotProduct(start, end, x_ncol, z_ncol, xtAy_local, A.isDotProductOptimized); // x'*A*y
 }
 MPI_Allreduce(&xtAy_local, &xtAy, 1, MPI_DOUBLE, MPI_SUM,
     MPI_COMM_WORLD);
 xtAy_local = 0.0;
 if (ierr) HPCG_fout << "Error in call to dot: " << ierr << ".\n" << endl;

 // Next, compute y'*A*x
 for (int i = 0; i < numBlocks; i++) {
   size_t start = blkInfo[i].start;
   size_t end = blkInfo[i].end;
   ComputeDotProduct(start, end, x_ncol, x_ncol, xNorm2_local, A.isDotProductOptimized);
 }
 MPI_Allreduce(&xNorm2_local, &xNorm2, 1, MPI_DOUBLE, MPI_SUM,
     MPI_COMM_WORLD);
 xNorm2_local = 0.0;
 #pragma omp parallel
 #pragma omp single
 ExchangeHalo(A,x_ncol,numNeighComms);
 #pragma omp taskwait
 for (int i = 0; i < numBlocks; i++) {
   size_t start = blkInfo[i].start;
   size_t end = blkInfo[i].end;
   ierr = ComputeSPMV(A, start, end, x_ncol, z_ncol); // b_computed = A*x_overlap
   if (ierr) HPCG_fout << "Error in call to SpMV: " << ierr << ".\n" << endl;
 }
 double ytAx = 0.0;
 double ytAx_local = 0.0;
 for (int i = 0; i < numBlocks; i++) {
   size_t start = blkInfo[i].start;
   size_t end = blkInfo[i].end;
   ierr = ComputeDotProduct(start, end, y_ncol, z_ncol, ytAx_local, A.isDotProductOptimized); // y'*A*x
 }
 MPI_Allreduce(&ytAx_local, &ytAx, 1, MPI_DOUBLE, MPI_SUM,
     MPI_COMM_WORLD);
 ytAx_local = 0.0;
 if (ierr) HPCG_fout << "Error in call to dot: " << ierr << ".\n" << endl;

 testsymmetry_data.depsym_spmv = std::fabs((long double) (xtAy - ytAx))/((xNorm2*ANorm*yNorm2 + yNorm2*ANorm*xNorm2) * (DBL_EPSILON));
 if (testsymmetry_data.depsym_spmv > 1.0) ++testsymmetry_data.count_fail;  // If the difference is > 1, count it wrong
 if (A.geom->rank==0) HPCG_fout << "Departure from symmetry (scaled) for SpMV abs(x'*A*y - y'*A*x) = " << testsymmetry_data.depsym_spmv << endl;

 // Test symmetry of multi-grid

 // Compute x'*Minv*y
 #pragma omp parallel
 #pragma omp single
 ierr = ComputeMG(A, y_ncol, z_ncol, numNeighComms); // z_ncol = Minv*y_ncol
 #pragma omp taskwait
 if (ierr) HPCG_fout << "Error in call to MG: " << ierr << ".\n" << endl;
 double xtMinvy = 0.0;
 double xtMinvy_local = 0.0;
 for (int i = 0; i < numBlocks; i++) {
   size_t start = blkInfo[i].start;
   size_t end = blkInfo[i].end;
   ierr = ComputeDotProduct(start, end, x_ncol, z_ncol, xtMinvy_local, A.isDotProductOptimized); // x'*Minv*y
 }
 MPI_Allreduce(&xtMinvy_local, &xtMinvy, 1, MPI_DOUBLE, MPI_SUM,
     MPI_COMM_WORLD);
 xtMinvy_local = 0.0;
 if (ierr) HPCG_fout << "Error in call to dot: " << ierr << ".\n" << endl;

 // Next, compute z'*Minv*x
 #pragma omp parallel
 #pragma omp single
 ierr = ComputeMG(A, x_ncol, z_ncol, numNeighComms); // z_ncol = Minv*x_ncol
 #pragma omp taskwait
 if (ierr) HPCG_fout << "Error in call to MG: " << ierr << ".\n" << endl;
 double ytMinvx = 0.0;
 double ytMinvx_local = 0.0;
 for (int i = 0; i < numBlocks; i++) {
   size_t start = blkInfo[i].start;
   size_t end = blkInfo[i].end;
   ierr = ComputeDotProduct(start, end, y_ncol, z_ncol, ytMinvx_local, A.isDotProductOptimized); // y'*Minv*x
 }
 MPI_Allreduce(&ytMinvx_local, &ytMinvx, 1, MPI_DOUBLE, MPI_SUM,
     MPI_COMM_WORLD);
 ytMinvx_local = 0.0;
 if (ierr) HPCG_fout << "Error in call to dot: " << ierr << ".\n" << endl;

 testsymmetry_data.depsym_mg = std::fabs((long double) (xtMinvy - ytMinvx))/((xNorm2*ANorm*yNorm2 + yNorm2*ANorm*xNorm2) * (DBL_EPSILON));
 if (testsymmetry_data.depsym_mg > 1.0) ++testsymmetry_data.count_fail;  // If the difference is > 1, count it wrong
 if (A.geom->rank==0) HPCG_fout << "Departure from symmetry (scaled) for MG abs(x'*Minv*y - y'*Minv*x) = " << testsymmetry_data.depsym_mg << endl;

 CopyVector(xexact, x_ncol); // Copy exact answer into overlap vector

 int numberOfCalls = 2;
 double residual = 0.0;
 #pragma omp parallel
 #pragma omp single
 ExchangeHalo(A,x_ncol,numNeighComms);
 #pragma omp taskwait
 for (int i=0; i< numberOfCalls; ++i) {
   for (int i = 0; i < numBlocks; i++) {
     size_t start = blkInfo[i].start;
     size_t end = blkInfo[i].end;
     ierr = ComputeSPMV(A, start, end, x_ncol, z_ncol); // b_computed = A*x_overlap
     if (ierr) HPCG_fout << "Error in call to SpMV: " << ierr << ".\n" << endl;
   }
   if ((ierr = ComputeResidual(A.localNumberOfRows, b, z_ncol, residual)))
     HPCG_fout << "Error in call to compute_residual: " << ierr << ".\n" << endl;
   if (A.geom->rank==0) HPCG_fout << "SpMV call [" << i << "] Residual [" << residual << "]" << endl;
 }
 DeleteVector(x_ncol);
 DeleteVector(y_ncol);
 DeleteVector(z_ncol);

 return 0;
}

