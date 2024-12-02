
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
 @file SetupHalo_ref.cpp

 HPCG routine
 */

#include <mpi.h>
#include <map>
#include <set>

#ifdef HPCG_DETAILED_DEBUG
#include <fstream>
using std::endl;
#include "hpcg.hpp"
#include <cassert>
#endif

#include "SetupHalo_ref.hpp"
#include "mytimer.hpp"

/*!
  Reference version of SetupHalo that prepares system matrix data structure and creates data necessary
  for communication of boundary values of this process.

  @param[inout] A    The known system matrix

  @see ExchangeHalo
*/
void SetupHalo_ref(SparseMatrix & A) {

  // Extract Matrix pieces

  local_int_t localNumberOfRows = A.localNumberOfRows;
  local_int_t ** mtxIndL = A.mtxIndL;

  // Scan global IDs of the nonzeros in the matrix.  Determine if the column ID matches a row ID.  If not:
  // 1) We call the ComputeRankOfMatrixRow function, which tells us the rank of the processor owning the row ID.
  //  We need to receive this value of the x vector during the halo exchange.
  // 2) We record our row ID since we know that the other processor will need this value from us, due to symmetry.

  std::map< int, std::set< std::pair<global_int_t, local_int_t> > > sendList;
  std::map< int, std::set< global_int_t> > receiveList;
  // std::map< int, std::set< global_int_t> > receiveListTmp;
  typedef std::map< int, std::set< std::pair<global_int_t, local_int_t> > >::iterator send_map_iter;
  typedef std::map< int, std::set< global_int_t> >::iterator recv_map_iter;
  typedef std::set< std::pair<global_int_t, local_int_t> >::iterator send_set_iter;
  typedef std::set<global_int_t>::iterator recv_set_iter;
  std::map< global_int_t, local_int_t > externalToLocalMap;

  global_int_t nx = A.geom->nx;
  global_int_t ny = A.geom->ny;
  global_int_t nz = A.geom->nz;
  global_int_t gnx = A.geom->gnx;
  global_int_t gny = A.geom->gny;
  global_int_t gnz = A.geom->gnz;
  global_int_t gix0 = A.geom->gix0;
  global_int_t giy0 = A.geom->giy0;
  global_int_t giz0 = A.geom->giz0;

  // TODO: maybe we can improve a bit with
  //       a parallel for critical here?
  for (local_int_t iz=0; iz<nz; iz++) {
    global_int_t giz = giz0+iz;
    for (local_int_t iy=0; iy<ny; iy++) {
      global_int_t giy = giy0+iy;
      for (local_int_t ix=0; ix<nx; ix++) {
        global_int_t gix = gix0+ix;
        if (!(iz == 0 || iz == nz -1 || iy == 0 || iy == ny -1 || ix == 0 || ix == nx -1))
          continue;
        local_int_t currentLocalRow = iz*nx*ny+iy*nx+ix;
        global_int_t currentGlobalRow = giz*gnx*gny+giy*gnx+gix;
        for (int sz=-1; sz<=1; sz++) {
          if (giz+sz>-1 && giz+sz<gnz) {
            for (int sy=-1; sy<=1; sy++) {
              if (giy+sy>-1 && giy+sy<gny) {
                for (int sx=-1; sx<=1; sx++) {
                  if (gix+sx>-1 && gix+sx<gnx) {
                    global_int_t curcol = currentGlobalRow+sz*gnx*gny+sy*gnx+sx;
#ifdef HPCG_DETAILED_DEBUG
                    local_int_t curlocalcol = currentLocalRow+sz*nx*ny+sy*nx+sx;
                    HPCG_fout << "rank, row , col, globalToLocalMap[col] = " << A.geom->rank << " " << currentGlobalRow << " "
                        << curcol << " " << curlocalcol << endl;
#endif
                    int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(A.geom), curcol);
                    if (A.geom->rank!=rankIdOfColumnEntry) {// If column index is not a row index, then it comes from another processor
                      receiveList[rankIdOfColumnEntry].insert(curcol);
                      sendList[rankIdOfColumnEntry].insert({currentGlobalRow, currentLocalRow}); // Matrix symmetry means we know the neighbor process wants my value
                    }
                  } // end x bounds test
                } // end sx loop
              } // end y bounds test
            } // end sy loop
          } // end z bounds test
        } // end sz loop
      } // end ix loop
    } // end iy loop
  } // end iz loop

  // Count number of matrix entries to send and receive
  local_int_t totalToBeSent = 0;
  for (send_map_iter curNeighbor = sendList.begin(); curNeighbor != sendList.end(); ++curNeighbor) {
    totalToBeSent += (curNeighbor->second).size();
  }
#ifdef HPCG_DETAILED_DEBUG
  local_int_t totalToBeReceived = 0;
  for (recv_map_iter curNeighbor = receiveList.begin(); curNeighbor != receiveList.end(); ++curNeighbor) {
    totalToBeReceived += (curNeighbor->second).size();
  }

  // These are all attributes that should be true, due to symmetry
  HPCG_fout << "totalToBeSent = " << totalToBeSent << " totalToBeReceived = " << totalToBeReceived << endl;
  assert(totalToBeSent==totalToBeReceived); // Number of sent entry should equal number of received
  assert(sendList.size()==receiveList.size()); // Number of send-to neighbors should equal number of receive-from
  // Each receive-from neighbor should be a send-to neighbor, and send the same number of entries
  for (recv_map_iter curNeighbor = receiveList.begin(); curNeighbor != receiveList.end(); ++curNeighbor) {
    assert(sendList.find(curNeighbor->first)!=sendList.end());
    assert(sendList[curNeighbor->first].size()==receiveList[curNeighbor->first].size());
  }
#endif

  // Build the arrays and lists needed by the ExchangeHalo function.
  double * sendBuffer = new double[totalToBeSent];
  local_int_t * elementsToSend = new local_int_t[totalToBeSent];
  int * neighbors = new int[sendList.size()];
  local_int_t * receiveLength = new local_int_t[receiveList.size()];
  local_int_t * sendLength = new local_int_t[sendList.size()];
  int neighborCount = 0;
  local_int_t receiveEntryCount = 0;
  local_int_t sendEntryCount = 0;
  for (recv_map_iter curNeighbor = receiveList.begin(); curNeighbor != receiveList.end(); ++curNeighbor, ++neighborCount) {
    int neighborId = curNeighbor->first; // rank of current neighbor we are processing
    neighbors[neighborCount] = neighborId; // store rank ID of current neighbor
    receiveLength[neighborCount] = receiveList[neighborId].size();
    sendLength[neighborCount] = sendList[neighborId].size(); // Get count if sends/receives
    for (recv_set_iter i = receiveList[neighborId].begin(); i != receiveList[neighborId].end(); ++i, ++receiveEntryCount) {
      externalToLocalMap[*i] = localNumberOfRows + receiveEntryCount; // The remote columns are indexed at end of internals
    }
    for (send_set_iter i = sendList[neighborId].begin(); i != sendList[neighborId].end(); ++i, ++sendEntryCount) {
      elementsToSend[sendEntryCount] = i->second; // store local ids of entry to send
    }
  }

/* Broken clang: https://github.com/llvm/llvm-project/issues/100536 */
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-compare"

  // Convert matrix indices to local IDs
  #pragma omp parallel
  #pragma omp single
  #pragma omp taskloop \
      default(shared)
  for (local_int_t iz=0; iz<nz; iz++) {
    global_int_t giz = giz0+iz;
    for (local_int_t iy=0; iy<ny; iy++) {
      global_int_t giy = giy0+iy;
      for (local_int_t ix=0; ix<nx; ix++) {
        global_int_t gix = gix0+ix;
        if (!(iz == 0 || iz == nz -1 || iy == 0 || iy == ny -1 || ix == 0 || ix == nx -1))
          continue;
        local_int_t currentLocalRow = iz*nx*ny+iy*nx+ix;
        global_int_t currentGlobalRow = giz*gnx*gny+giy*gnx+gix;
        local_int_t * currentIndexPointerL = mtxIndL[currentLocalRow]; // Pointer to current index in current row
        for (int sz=-1; sz<=1; sz++) {
          if (giz+sz>-1 && giz+sz<gnz) {
            for (int sy=-1; sy<=1; sy++) {
              if (giy+sy>-1 && giy+sy<gny) {
                for (int sx=-1; sx<=1; sx++) {
                  if (gix+sx>-1 && gix+sx<gnx) {
                    global_int_t curcol = currentGlobalRow+sz*gnx*gny+sy*gnx+sx;
                    int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(A.geom), curcol);
                    // My column index are already set up
                    if (A.geom->rank != rankIdOfColumnEntry) {
                      // If column index is not a row index, then it comes from another processor
                      *currentIndexPointerL = externalToLocalMap[curcol];
                    }
                    currentIndexPointerL++;
                  } // end x bounds test
                } // end sx loop
              } // end y bounds test
            } // end sy loop
          } // end z bounds test
        } // end sz loop
      } // end ix loop
    } // end iy loop
  } // end iz loop

#pragma clang diagnostic pop

  // Store contents in our matrix struct
  A.numberOfExternalValues = externalToLocalMap.size();
  A.localNumberOfColumns = A.localNumberOfRows + A.numberOfExternalValues;
  A.numberOfSendNeighbors = sendList.size();
  A.totalToBeSent = totalToBeSent;
  A.elementsToSend = elementsToSend;
  A.neighbors = neighbors;
  A.receiveLength = receiveLength;
  A.sendLength = sendLength;
  A.sendBuffer = sendBuffer;

#ifdef HPCG_DETAILED_DEBUG
  HPCG_fout << " For rank " << A.geom->rank << " of " << A.geom->size << ", number of neighbors = " << A.numberOfSendNeighbors << endl;
  for (int i = 0; i < A.numberOfSendNeighbors; i++) {
    HPCG_fout << "     rank " << A.geom->rank << " neighbor " << neighbors[i] << " send/recv length = " << sendLength[i] << "/" << receiveLength[i] << endl;
    for (local_int_t j = 0; j<sendLength[i]; ++j)
      HPCG_fout << "       rank " << A.geom->rank << " elementsToSend[" << j << "] = " << elementsToSend[j] << endl;
  }
#endif

  return;
}
