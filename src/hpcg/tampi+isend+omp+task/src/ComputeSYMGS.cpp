
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
 @file ComputeSYMGS.cpp

 HPCG routine
 */

#include "ComputeSYMGS.hpp"
#include "ComputeSYMGS_ref.hpp"

/*!
  Routine to compute one step of symmetric Gauss-Seidel:

  Assumption about the structure of matrix A:
  - Each row 'i' of the matrix has nonzero diagonal value whose address is matrixDiagonal[i]
  - Entries in row 'i' are ordered such that:
       - lower triangular terms are stored before the diagonal element.
       - upper triangular terms are stored after the diagonal element.
       - No other assumptions are made about entry ordering.

  Symmetric Gauss-Seidel notes:
  - We use the input vector x as the RHS and start with an initial guess for y of all zeros.
  - We perform one forward sweep.  Since y is initially zero we can ignore the upper triangular terms of A.
  - We then perform one back sweep.
       - For simplicity we include the diagonal contribution in the for-j loop, then correct the sum after

  @param[in] A the known system matrix
  @param[in] r the input vector
  @param[inout] x On entry, x should contain relevant values, on exit x contains the result of one symmetric GS sweep with r as the RHS.

  @return returns 0 upon success and non-zero otherwise

  @warning Early versions of this kernel (Version 1.1 and earlier) had the r and x arguments in reverse order, and out of sync with other kernels.

  @see ComputeSYMGS_ref
*/

int ComputeSYMGS( const SparseMatrix & A, const Vector & r, Vector & x) {
  double *const &rv = r.values;
  double *const &xv = x.values;

  const BlockLayout &blkLayout = A.blkLayout;
  const std::vector< std::vector<local_int_t> > &colorBlkInfo = blkLayout.colorBlkInfo;
  const std::vector<BlockInfo> &blkInfo = blkLayout.blkInfo;
  const std::vector< std::vector<local_int_t> > &blkDepList = blkLayout.blkDepList;

  for (int k = 0; k < (int) colorBlkInfo.size(); ++k) { // Color loop
    for (int j = 0; j < (int) colorBlkInfo[k].size(); j++) { // Block loop
      int block = colorBlkInfo[k][j];
      int start = blkInfo[block].start;
      int end = blkInfo[block].end;
      #pragma omp task \
          default(shared) \
          firstprivate(start, end) \
          depend(iterator(d=0 : blkDepList[block].size()), in: \
                 rv[blkInfo[blkDepList[block][d]].start] ) \
          depend(iterator(d=0 : blkDepList[block].size()), in: \
                 xv[blkInfo[blkDepList[block][d]].start] ) \
          depend(out: xv[start] )
      {
        // in(A, r) inout(x)
        ComputeFGS(A, start, end, r, x);
      }
    }
  }
  for (int k = colorBlkInfo.size() - 1; k >= 0; --k) { // Color loop
    for (int j = 0; j < (int) colorBlkInfo[k].size(); j++) { // Block loop
      int block = colorBlkInfo[k][j];
      int start = blkInfo[block].start;
      int end = blkInfo[block].end;
      #pragma omp task \
          default(shared) \
          firstprivate(start, end) \
          depend(iterator(d=0 : blkDepList[block].size()), in: \
                 rv[blkInfo[blkDepList[block][d]].start] ) \
          depend(iterator(d=0 : blkDepList[block].size()), in: \
                 xv[blkInfo[blkDepList[block][d]].start] ) \
          depend(out: xv[start] )
      {
        // in(A, r) inout(x)
        ComputeBGS(A, start, end, r, x);
      }
    }
  }
  return 0;
}

int ComputeFGS( const SparseMatrix & A, local_int_t start, local_int_t end, const Vector & r, Vector & x) {

  assert(x.localLength==A.localNumberOfColumns); // Make sure x contain space for halo values

  double ** matrixDiagonal = A.matrixDiagonal;  // An array of pointers to the diagonal entries A.matrixValues
  const double * const rv = r.values;
  double * const xv = x.values;

  for (local_int_t i=start; i< end; i++) {
    const double * const currentValues = A.matrixValues[i];
    const local_int_t * const currentColIndices = A.mtxIndL[i];
    const int currentNumberOfNonzeros = A.nonzerosInRow[i];
    const double  currentDiagonal = matrixDiagonal[i][0]; // Current diagonal value
    double sum = rv[i]; // RHS value

    for (int j=0; j< currentNumberOfNonzeros; j++) {
      local_int_t curCol = currentColIndices[j];
      sum -= currentValues[j] * xv[curCol];
    }
    sum += xv[i]*currentDiagonal; // Remove diagonal contribution from previous loop

    xv[i] = sum/currentDiagonal;

  }

  return 0;

}

int ComputeBGS( const SparseMatrix & A, local_int_t start, local_int_t end, const Vector & r, Vector & x) {

  assert(x.localLength==A.localNumberOfColumns); // Make sure x contain space for halo values
  assert(end <= A.localNumberOfColumns); // Make sure x contain space for halo values

  double ** matrixDiagonal = A.matrixDiagonal;  // An array of pointers to the diagonal entries A.matrixValues
  const double * const rv = r.values;
  double * const xv = x.values;

  // Now the back sweep.

  for (local_int_t i=end-1; i>=start; i--) {
    const double * const currentValues = A.matrixValues[i];
    const local_int_t * const currentColIndices = A.mtxIndL[i];
    const int currentNumberOfNonzeros = A.nonzerosInRow[i];
    const double  currentDiagonal = matrixDiagonal[i][0]; // Current diagonal value
    double sum = rv[i]; // RHS value

    for (int j = 0; j< currentNumberOfNonzeros; j++) {
      local_int_t curCol = currentColIndices[j];
      sum -= currentValues[j]*xv[curCol];
    }
    sum += xv[i]*currentDiagonal; // Remove diagonal contribution from previous loop

    xv[i] = sum/currentDiagonal;
  }

  return 0;

}
