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

#include <assert.h>
#include <mpi.h>

#include "block.h"
#include "comm.h"
#include "proto.h"

// Pack and unpack a block to move blocks between processors.
void pack_block(int n, num_sz number, int level, double *buffer)
{
   // Pack a block's values into a buffer to send. This and the
   // next routine which unpacks the values depend on the assumption
   // that the size of an integer is no larger than a double.
   // Note: number and level fields of the block should not be
   // accessed since they may be modified by the exchange function
   // before calling this pack function. Instead, copy the number
   // and level parameters into the send buffer

   int v, i, j, k, l;
   int *send_int = (int *) buffer;
   long long *send_ll = (long long *) buffer;
   block *bp = &blocks[n];

   send_ll[0] = (long long) number;
   if (bp->parent_node == my_pe && bp->parent != -1)
      // bp->parent converted to parent number (from index) in move_blocks
      send_ll[1] = (long long) (-2 - bp->parent);
   else
      send_ll[1] = (long long) bp->parent;
   l = 4;
   send_int[l++] = level;
   send_int[l++] = bp->refine;
   send_int[l++] = bp->parent_node;
   send_int[l++] = bp->child_number;
   for (i = 0; i < 6; i++) {
      send_int[l++] = bp->nei_refine[i];
      send_int[l++] = bp->nei_level[i];
      for (j = 0; j < 2; j++)
         for (k = 0; k < 2; k++)
            send_int[l++] = bp->nei[i][j][k];
   }
   for (i = 0; i < 3; i++)
      send_int[l++] = bp->cen[i];

   typedef double (*block3d_t)[y_block_size+2][z_block_size+2];

   for (v = 0; v < num_vars; v++) {
      block3d_t matrix = (block3d_t) &(bp->array[v*block3d_size]);
      for (i = 1; i <= x_block_size; i++)
         for (j = 1; j <= y_block_size; j++)
            for (k = 1; k <= z_block_size; k++)
               buffer[l++] = matrix[i][j][k];
   }
}

void unpack_block(int n, double *buffer)
{
   int v, i, j, k, l;
   int *recv_int = (int *) buffer;
   long long *recv_ll = (long long *) buffer;
   block *bp = &blocks[n];

   // These fields should not be modified since
   // were already set by the exchange function
   assert(bp->new_proc == -1);
   assert(bp->number == (num_sz) recv_ll[0]);
   assert(bp->level == recv_int[4]);

   bp->parent = (num_sz) recv_ll[1];
   l = 5;
   bp->refine = recv_int[l++];
   bp->parent_node = recv_int[l++];
   bp->child_number = recv_int[l++];
   for (i = 0; i < 6; i++) {
      bp->nei_refine[i] = recv_int[l++];
      bp->nei_level[i] = recv_int[l++];
      for (j = 0; j < 2; j++)
         for (k = 0; k < 2; k++)
            bp->nei[i][j][k] = recv_int[l++];
   }
   for (i = 0; i < 3; i++)
      bp->cen[i] = recv_int[l++];

   typedef double (*block3d_t)[y_block_size+2][z_block_size+2];

   for (v = 0; v < num_vars; v++) {
      block3d_t matrix = (block3d_t) &(bp->array[v*block3d_size]);
      for (i = 1; i <= x_block_size; i++)
         for (j = 1; j <= y_block_size; j++)
            for (k = 1; k <= z_block_size; k++)
               matrix[i][j][k] = buffer[l++];
   }
}
