
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
 @file ExchangeHalo.cpp

 HPCG routine
 */

// Compile this routine only if running with MPI
#include <mpi.h>
#include "Geometry.hpp"
#include "ExchangeHalo.hpp"
#include <cstdlib>
#include <iostream>
#include <TAMPI.h>

/*!
  Communicates data that is at the border of the part of the domain assigned to this processor.

  @param[in]    A The known system matrix
  @param[inout] x On entry: the local vector entries followed by entries to be communicated; on exit: the vector with non-local entries updated by other processors
 */

void ExchangeHalo(const SparseMatrix & A, Vector & x, int numNeighComms) {

  // Extract Matrix pieces

  local_int_t localNumberOfRows = A.localNumberOfRows;
  int * neighbors = A.neighbors;
  double * sendBuffer = A.sendBuffer;
  local_int_t * elementsToSend = A.elementsToSend;

  double * const xv = x.values;

  int size, rank; // Number of MPI processes, My process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int MPI_MY_TAG = 9999;
  // TODO: assign unique MPI_TAG. Use a cleaner way to do this...
  static int n_exchange = 0;

  //
  // Externals are at end of locals
  //
  double * x_external = (double *) xv + localNumberOfRows;

  //
  // Fill up send buffer
  //

  const BlockLayout & blkLayout = A.blkLayout;
  const std::vector< std::map< local_int_t, std::vector<local_int_t> > > &neighborHaloToBlkHalo = blkLayout.neighborHaloToBlkHalo;

  for (size_t k = 0; k < neighborHaloToBlkHalo.size(); ++k) {
    std::map<local_int_t, std::vector<local_int_t> >::const_iterator it;
    int extraBlockIdx = 0;
    for (it = neighborHaloToBlkHalo[k].begin();
         it != neighborHaloToBlkHalo[k].end(); ++it) {

      assert(extraBlockIdx < numNeighComms);

      std::vector<local_int_t> const& list_idx = it->second;
      int size = list_idx.size();
      #pragma omp task shared(list_idx) firstprivate(n_exchange) \
                       depend(iterator(elem=0:size), in: xv[ list_idx[elem] ] ) \
                       depend(out: sendBuffer[0] )
      {
        for (local_int_t i=0; i<size; i++) sendBuffer[i] = xv[elementsToSend[i]];
#if (TAMPI_VERSION_MAJOR < 4)
        MPI_Request request;
        MPI_Isend(sendBuffer, size, MPI_DOUBLE, neighbors[k], MPI_MY_TAG + n_exchange*numNeighComms + extraBlockIdx, MPI_COMM_WORLD, &request);
        TAMPI_Iwaitall(1, &request, MPI_STATUSES_IGNORE);
#else
        TAMPI_Isend(sendBuffer, size, MPI_DOUBLE, neighbors[k], MPI_MY_TAG + n_exchange*numNeighComms + extraBlockIdx, MPI_COMM_WORLD);
#endif
      }

      ++extraBlockIdx;
      sendBuffer += size;
      elementsToSend += size;
    }
  }

  for (size_t k = 0; k < neighborHaloToBlkHalo.size(); ++k) {
    std::map<local_int_t, std::vector<local_int_t> >::const_iterator it;

    int extraBlockIdx = 0;
    for (it = neighborHaloToBlkHalo[k].begin();
         it != neighborHaloToBlkHalo[k].end(); ++it) {

      assert(extraBlockIdx < numNeighComms);

      std::vector<local_int_t> const& list_idx = it->second;
      int size = list_idx.size();

      #pragma omp task firstprivate(n_exchange) depend(out: x_external[0] )
      {
#if (TAMPI_VERSION_MAJOR < 4)
        MPI_Request request;
        MPI_Irecv(x_external, size, MPI_DOUBLE, neighbors[k], MPI_MY_TAG + n_exchange*numNeighComms + extraBlockIdx, MPI_COMM_WORLD, &request);
        TAMPI_Iwaitall(1, &request, MPI_STATUSES_IGNORE);
#else
        TAMPI_Irecv(x_external, size, MPI_DOUBLE, neighbors[k], MPI_MY_TAG + n_exchange*numNeighComms + extraBlockIdx, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
      }
      x_external += size;
      ++extraBlockIdx;
    }
  }
  ++n_exchange;

  return;
}
