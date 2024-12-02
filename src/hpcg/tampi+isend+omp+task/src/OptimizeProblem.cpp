
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
 @file OptimizeProblem.cpp

 HPCG routine
 */

#include <mpi.h>

#include "hpcg.hpp"
#include "OptimizeProblem.hpp"

#ifdef HPCG_DEBUG
#include <iostream>
static void PrintBlockInfo(SparseMatrix &A) {
  const BlockLayout &blkLayout = A.blkLayout;
  local_int_t numBlocks = blkLayout.numBlocks;
  const std::vector<BlockInfo> &blkInfo = blkLayout.blkInfo;
  const std::vector< std::vector<local_int_t> > &blkDepList = blkLayout.blkDepList;
  const std::vector< std::vector<local_int_t> > &colorBlkInfo = blkLayout.colorBlkInfo;
  const std::vector< std::vector<local_int_t> > &coarserToCurrentBlkInfo = blkLayout.coarserToCurrentBlkInfo;
  for (int i = 0; i < numBlocks; ++i) {
    // Print Block ranges
    std::cout
    << "BLOCK " << i
    << " [" << blkInfo[i].start << ", " << blkInfo[i].end
    << "]:";
    const local_int_t * const cur_inds = A.mtxIndL[i];
    // and Block dependencies
    for (int g1 = 0; g1 < blkDepList[i].size(); g1++) {
      std::cout << " " << blkDepList[i][g1];
    }
    std::cout << "\n";
  }
  for (int i = numBlocks; i < blkInfo.size(); ++i) {
    std::cout
    << "BLOCK " << i
    << " [" << blkInfo[i].start << ", " << blkInfo[i].end
    << "]\n";
  }
  for (int i = 0; i < colorBlkInfo.size(); ++i) {
      std::cout << "COLOR " << i << ":";
      for (int j = 0; j < colorBlkInfo[i].size(); ++j) {
          std::cout << " " << colorBlkInfo[i][j];
      }
      std::cout << std::endl;
  }

  for (int i = 0; i < coarserToCurrentBlkInfo.size(); ++i) {
    std::cout << "CoarseLevel " << i << " uses :";
    for (int j = 0; j < coarserToCurrentBlkInfo[i].size(); ++j) {
      std::cout << " " << coarserToCurrentBlkInfo[i][j];
    }
    std::cout << std::endl;
  }
  // TODO: print neightbors info
}
#endif

/*!
  Optimizes the data structures used for CG iteration to increase the
  performance of the benchmark version of the preconditioned CG algorithm.

  @param[inout] A      The known system matrix, also contains the MG hierarchy in attributes Ac and mgData.
  @param[inout] data   The data structure with all necessary CG vectors preallocated
  @param[inout] b      The known right hand side vector
  @param[inout] x      The solution vector to be computed in future CG iteration
  @param[inout] xexact The exact solution vector

  @return returns 0 upon success and non-zero otherwise

  @see GenerateGeometry
  @see GenerateProblem
*/
int OptimizeProblem(
    SparseMatrix & A, CGData & data, Vector & b, Vector & x, Vector & xexact,
    int numBlocks, int numNeighComms) {

  // This function can be used to completely transform any part of the data structures.
  // Right now it does nothing, so compiling with a check for unused variables results in complaints

  // Fill block info

  die_if(numBlocks <= 0, "Expected --nblocks=<num_blocks> greater than 0");
  die_if(numNeighComms <= 0, "Expected --ncomms=<num_neigh_comms> greater than 0");
  int chunk = A.localNumberOfRows/(int)(numBlocks);
  // OptimizeProblem is called one time for each multigrid level.
  static int level = 0;
  if (A.geom->rank==0)
    printf("blocksize %d %lu\n", level, chunk*sizeof(double));
  level++;
  die_if(!chunk, "chunk == 0");
  die_if(numNeighComms > 1 && A.localNumberOfRows%(A.geom->nx*A.geom->ny), "Using ncomms > 1 expects chunk to be multiple of nx*ny");

  std::vector<BlockInfo> blkInfo(numBlocks);
  for (int i = 0; i < numBlocks; ++i) {
    blkInfo[i].start = i * chunk;
    blkInfo[i].end = (i + 1) * chunk;
    if (i == numBlocks - 1) {
      blkInfo[i].end = A.localNumberOfRows;
    }
  }

  int numNeighbors = A.numberOfSendNeighbors;
  local_int_t * receiveLength = A.receiveLength;
  local_int_t * elementsToSend = A.elementsToSend;

  std::vector< std::map< local_int_t, std::vector<local_int_t> > > neighborHaloToBlkHalo(numNeighbors);

  local_int_t acc = 0;
  local_int_t extraBlocks = 0;
  for (int nbor = 0; nbor < numNeighbors; nbor++) {
    local_int_t n_recv = receiveLength[nbor];
    local_int_t comm_chunk = n_recv/numNeighComms;
    if (comm_chunk == 0) comm_chunk = 1;
    int n_block = 0;
    for (int ircv = 0; ircv < n_recv; ircv += comm_chunk) {
      local_int_t start = ircv;
      local_int_t end = start + comm_chunk > n_recv ? n_recv : start + comm_chunk;
      for (int iircv = start; iircv < end; ++iircv) {
        // Normalize element to send idx to block
        // TODO: Right now we keep repeated indexes because
        // the size of the array is used in ExchangeHalo. We may want
        // to materialize how many indices we have and remove repeateds
        int iircv_norm = (elementsToSend[acc + iircv]/chunk)*chunk;
        neighborHaloToBlkHalo[nbor][n_block].push_back(iircv_norm);
      }
      assert(start != end);
      ++n_block;
      blkInfo.push_back({A.localNumberOfRows + acc + start,
                         A.localNumberOfRows + acc + end});
      ++extraBlocks;
    }
    acc += n_recv;
  }

  std::vector< std::vector<local_int_t> > blkDepList;
  std::vector<bool> accessed_inds(A.localNumberOfColumns, false);
  for (int i = 0; i < (int) numBlocks; i++) {
    int start = blkInfo[i].start;
    int end = blkInfo[i].end;
    for (int j = start; j < end; j++) {
      const local_int_t * const cur_inds = A.mtxIndL[j];
      const int cur_nnz = A.nonzerosInRow[j];

      for (int k = 0; k < cur_nnz; k++)
        accessed_inds[cur_inds[k]] = true;
    }
    for (int j1 = 0; j1 < numBlocks + extraBlocks; ++j1) {
      int start1 = blkInfo[j1].start;
      int end1 = blkInfo[j1].end;
      for (int k = start1; k < end1; k++) {
        if (accessed_inds[k]) {
          if (i >= (int) blkDepList.size())
            blkDepList.push_back(std::vector<local_int_t>());
          blkDepList[i].push_back(j1);
          break;
        }
      }
    }
    std::fill(accessed_inds.begin(), accessed_inds.end(), false);
  }

  // Make a table of assignable colors for each block
  std::vector< std::vector<bool> > possibleColors(
      numBlocks, std::vector<bool> (numBlocks, true));
  std::vector< std::vector<local_int_t> > colorBlkInfo;
  std::vector<int> color_freq(numBlocks, 0);
  const int MAX_BLOCKS_PER_COLOR = -1;

  // Assign a color from among the possible ones and mark as used to the adjacent blocks.
  for (int i = 0; i < numBlocks; ++i) { // Block loop
    for (int j = 0; j < numBlocks; ++j) { // Color loop
      if (possibleColors[i][j]) {
        // Assign color
        if (j >= (int) colorBlkInfo.size()) {
          colorBlkInfo.push_back(std::vector<local_int_t>());
        }
        colorBlkInfo[j].push_back(i);
        // Mark all forward adjacents skipping
        // extraBlocks, or Mark all if we reached MAX_BLOCKS_PER_COLOR
        if (++color_freq[j] == MAX_BLOCKS_PER_COLOR) {
          for (int k = i; k < numBlocks; ++k)
            possibleColors[k][j] = false;
        } else {
          for (int k = 0; k < (int) blkDepList[i].size(); ++k) {
            if (blkDepList[i][k] > i
                && blkDepList[i][k] < numBlocks)
              possibleColors[blkDepList[i][k]][j] = false;
          }
        }
        break;
      }
    }
  }

  // If there are coarser levels, colorize them also
  if (A.mgData) {
    OptimizeProblem(*A.Ac, data, b, x, xexact, numBlocks, numNeighComms);

    const BlockLayout & blkCLayout = A.Ac->blkLayout;
    int numCBlocks = blkCLayout.numBlocks;
    const std::vector<BlockInfo> &blkCInfo = blkCLayout.blkInfo;

    local_int_t *f2c = A.mgData->f2cOperator;

    std::vector< std::vector<local_int_t> > coarserToCurrentBlkInfo(numCBlocks);
    for (int i = 0; i < numCBlocks; i++) {
      int start = blkCInfo[i].start;
      int end = blkCInfo[i].end;
      std::vector<bool> curBlocksUsed(numBlocks, false);
      for (local_int_t j = start; j < end; ++j) {
        int coarseInd = f2c[j];
        for (int k = 0; k < numBlocks; ++k) {
          int end1 = blkInfo[k].end;
          if (coarseInd < end1) {
            curBlocksUsed[k] = true;
            break;
          }
        }
      }
      for (int j = 0; j < (int) curBlocksUsed.size(); ++j) {
        if (curBlocksUsed[j])
          coarserToCurrentBlkInfo[i].push_back(j);
      }
    }
    A.blkLayout.coarserToCurrentBlkInfo = coarserToCurrentBlkInfo;
  }


  A.blkLayout.blkInfo = blkInfo;
  A.blkLayout.numBlocks = numBlocks;
  A.blkLayout.blkDepList = blkDepList;
  A.blkLayout.colorBlkInfo = colorBlkInfo;
  A.blkLayout.extraBlocks = extraBlocks;
  A.blkLayout.neighborHaloToBlkHalo = neighborHaloToBlkHalo;

#ifdef HPCG_DEBUG
  PrintBlockInfo(A);
#endif

  return 0;
}

// Helper function (see OptimizeProblem.hpp for details)
double OptimizeProblemMemoryUse(const SparseMatrix & A) {
  (void) A;

  return 0.0;

}
