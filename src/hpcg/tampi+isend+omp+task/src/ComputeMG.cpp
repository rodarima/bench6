
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
 @file ComputeMG.cpp

 HPCG routine
 */

#include "ComputeMG.hpp"
#include "ComputeMG_ref.hpp"

#include "ComputeSYMGS.hpp"
#include "ComputeSPMV.hpp"
#include "ComputeRestriction.hpp"
#include "ComputeProlongation.hpp"

#include "ExchangeHalo.hpp"

#include <iostream>
#include <mpi.h>

/*!
  @param[in] A the known system matrix
  @param[in] r the input vector
  @param[inout] x On exit contains the result of the multigrid V-cycle with r as the RHS, x is the approximation to Ax = r.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeMG_ref
*/

int ComputeMG(const SparseMatrix  & A, const Vector & r, Vector & x, int numNeighComms) {
  assert(x.localLength==A.localNumberOfColumns); // Make sure x contain space for halo values

  double *const &xv = x.values;

  const BlockLayout &blkLayout = A.blkLayout;
  const local_int_t &numBlocks = blkLayout.numBlocks;
  const std::vector<BlockInfo> &blkInfo = blkLayout.blkInfo;
  const std::vector< std::vector<local_int_t> > &blkDepList = blkLayout.blkDepList;

  for (int i = 0; i < numBlocks; i++) {
    int start = blkInfo[i].start;
    int end = blkInfo[i].end;
    #pragma omp task \
        default(shared) \
        firstprivate(i, start, end) \
        depend(out:xv[start] )
    ZeroVector(start, end, x); // initialize x to zero
  }

  if (A.mgData!=0) { // Go to next coarse level if defined
    const BlockLayout &blkCLayout = A.Ac->blkLayout;
    int numCBlocks = blkCLayout.numBlocks;
    const std::vector<BlockInfo> &blkCInfo = blkCLayout.blkInfo;
    const std::vector< std::vector<local_int_t> > &coarserToCurrentBlkInfo = blkLayout.coarserToCurrentBlkInfo;

    double *const &Axfv = A.mgData->Axf->values;
    double *const &rcv = A.mgData->rc->values;
    double *const &xcv = A.mgData->xc->values;

    int numberOfPresmootherSteps = A.mgData->numberOfPresmootherSteps;
    for (int i=0; i< numberOfPresmootherSteps; ++i) {
      ExchangeHalo(A,x,numNeighComms);
      ComputeSYMGS(A, r, x);
    }
    ExchangeHalo(A,x,numNeighComms);
    for (int i = 0; i < numBlocks; i++) {
      int start = blkInfo[i].start;
      int end = blkInfo[i].end;
      // in(A, x) out(*A.mgData->Axf)
      #pragma omp task \
          default(shared) \
          shared(blkInfo, blkDepList) \
          firstprivate(start, end) \
          depend(iterator(j=0:blkDepList[i].size()), in: xv[blkInfo[blkDepList[i][j]].start] ) \
          depend(out: Axfv[start] )
      ComputeSPMV(A, start, end, x, *A.mgData->Axf);
    }
    // Perform restriction operation using simple injection

    // in(r, *A.mgData->Axf) out(A.mgData->rc)
    for (int i = 0; i < numCBlocks; i++) {
      int start = blkCInfo[i].start;
      int end = blkCInfo[i].end;
      #pragma omp task \
          default(shared) \
          shared(blkInfo, coarserToCurrentBlkInfo) \
          firstprivate(start, end) \
          depend(iterator(j=0:coarserToCurrentBlkInfo[i].size()), in: Axfv[blkInfo[coarserToCurrentBlkInfo[i][j]].start] ) \
          depend(out: rcv[start] )
      ComputeRestriction(A, start, end, r);
    }

    // in(*A.Ac, *A.mgData->rc) inout(*A.mgData->xc)
    ComputeMG(*A.Ac,*A.mgData->rc, *A.mgData->xc, numNeighComms);

    // in(Af.mgData->xc) inout(x)
    for (int i = 0; i < numCBlocks; i++) {
      int start = blkCInfo[i].start;
      int end = blkCInfo[i].end;
        // depend(iterator(j=0:coarserToCurrentBlkInfo[i].size()), inoutset: xv[blkInfo[coarserToCurrentBlkInfo[i][j]].start] )
      #pragma omp task \
          default(shared) \
          shared(blkInfo, coarserToCurrentBlkInfo) \
          firstprivate(start, end) \
          depend(iterator(j=0:coarserToCurrentBlkInfo[i].size()), inout: xv[blkInfo[coarserToCurrentBlkInfo[i][j]].start] ) \
          depend(in: xcv[start] )
      ComputeProlongation(A, start, end, x);
    }

    int numberOfPostsmootherSteps = A.mgData->numberOfPostsmootherSteps;
    for (int i=0; i< numberOfPostsmootherSteps; ++i) {
      ExchangeHalo(A,x,numNeighComms);
      ComputeSYMGS(A, r, x);
    }
  }
  else {
    ExchangeHalo(A,x,numNeighComms);
    ComputeSYMGS(A, r, x);
  }
  return 0;
}
