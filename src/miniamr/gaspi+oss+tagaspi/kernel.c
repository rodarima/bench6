// ************************************************************************
//
// miniAMR: stencil computations with boundary exchange and AMR.
//
// Copyright (2014) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
// Questions? Contact Courtenay T. Vaughan (ctvaugh@sandia.gov)
//                    Richard F. Barrett (rfbarre@sandia.gov)
//
// ************************************************************************

//
// This file is part of miniAMR and is licensed under the terms contained in the ompss-2-tagaspi/LICENSE file.
//
// Copyright (C) 2020-2021 Barcelona Supercomputing Center (BSC)
//

#include <string.h>

#include "block.h"
#include "comm.h"
#include "proto.h"

void stencil_kernel(double *matrix_array, double *work_array)
{
   typedef double (*block3d_t)[y_block_size+2][z_block_size+2];
   block3d_t __restrict__ matrix = (block3d_t) matrix_array;
   block3d_t __restrict__ work = (block3d_t) work_array;

   size_t msize = (x_block_size+2)*(y_block_size+2)*(z_block_size+2);
   memcpy(work, matrix, msize*sizeof(double));

   for (int i = 1; i <= x_block_size; i++)
      for (int j = 1; j <= y_block_size; j++)
         for (int k = 1; k <= z_block_size; k++)
            matrix[i][j][k] = (work[i-1][j  ][k  ] +
                               work[i  ][j-1][k  ] +
                               work[i  ][j  ][k-1] +
                               work[i  ][j  ][k  ] +
                               work[i  ][j  ][k+1] +
                               work[i  ][j+1][k  ] +
                               work[i+1][j  ][k  ])/7.0;
}

void copy_ghosts_westeast(block *bp, block *bp1, int start, int num_comm)
{
   typedef double (*block3d_t)[y_block_size+2][z_block_size+2];
   for (int m = start; m < start+num_comm; m++) {
      block3d_t __restrict__ matrix = (block3d_t) &(bp->array[m*block3d_size]);
      block3d_t __restrict__ matrix1 = (block3d_t) &(bp1->array[m*block3d_size]);
      for (int j = 1; j <= y_block_size; j++)
         for (int k = 1; k <= z_block_size; k++) {
            matrix1[x_block_size+1][j][k] = matrix[1][j][k];
            matrix[0][j][k] = matrix1[x_block_size][j][k];
         }
   }
}

void copy_ghosts_southnorth(block *bp, block *bp1, int start, int num_comm, int is, int ie)
{
   typedef double (*block3d_t)[y_block_size+2][z_block_size+2];
   for (int m = start; m < start+num_comm; m++) {
      block3d_t __restrict__ matrix = (block3d_t) &(bp->array[m*block3d_size]);
      block3d_t __restrict__ matrix1 = (block3d_t) &(bp1->array[m*block3d_size]);
      for (int i = is; i <= ie; i++)
         for (int k = 1; k <= z_block_size; k++) {
            matrix1[i][y_block_size+1][k] = matrix[i][1][k];
            matrix[i][0][k] = matrix1[i][y_block_size][k];
         }
   }
}

void copy_ghosts_downup(block *bp, block *bp1, int start, int num_comm, int is, int ie, int js, int je)
{
   typedef double (*block3d_t)[y_block_size+2][z_block_size+2];
   for (int m = start; m < start+num_comm; m++) {
      block3d_t __restrict__ matrix = (block3d_t) &(bp->array[m*block3d_size]);
      block3d_t __restrict__ matrix1 = (block3d_t) &(bp1->array[m*block3d_size]);
      for (int i = is; i <= ie; i++)
         for (int j = js; j <= je; j++) {
            matrix1[i][j][z_block_size+1] = matrix[i][j][1];
            matrix[i][j][0] = matrix1[i][j][z_block_size];
         }
   }
}
