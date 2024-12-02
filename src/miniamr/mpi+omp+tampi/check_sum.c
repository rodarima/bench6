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
// This file is part of miniAMR and is licensed under the terms contained in the openmp-tasks/LICENSE file.
//
// Copyright (C) 2020-2021 Barcelona Supercomputing Center (BSC)
//

#include <math.h>
#include <mpi.h>
#include <stdio.h>

#include "block.h"
#include "comm.h"
#include "timer.h"
#include "proto.h"

// Generate check sum for a variable over all active blocks.
void check_sum(int var, int num, double *sum)
{
   check_sum_local(var, num, sum);
   //#pragma omp taskwait

   check_sum_remote(var, num, sum);

   total_red++;
}

// Generate local check sum for a variable over all active blocks.
void check_sum_local(int var, int num, double *sum)
{
   double t1 = timer();

   for (int in = 0; in < sorted_index[num_refine+1]; in++) {
      double *barray = blocks[sorted_list[in].n].array;

      //#pragma omp task depend(in: barray[var*block3d_size:num*block3d_size])
      for (int v = var; v < var+num; v++) {
         typedef double (*block3d_t)[y_block_size+2][z_block_size+2];
         block3d_t matrix = (block3d_t) &(barray[v*block3d_size]);

         double block_sum = 0.0;
         for (int i = 1; i <= x_block_size; i++)
            for (int j = 1; j <= y_block_size; j++)
               for (int k = 1; k <= z_block_size; k++)
                  block_sum += matrix[i][j][k];

         //#pragma omp atomic
         sum[v] += block_sum;
      }
   }
   timer_cs_calc += timer() - t1;
   total_red++;
}

// Reduce remotely the check sum for a variable over all active blocks.
void check_sum_remote(int var, int num, double *sum)
{
   double t1 = timer();
   MPI_Allreduce(MPI_IN_PLACE, &sum[var], num, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   timer_cs_red += timer() - t1;
}

void check_sum_verify(int var, int num, double *sum, double *grid_sum, int ts)
{
   for (int v = var; v < num; v++) {
      if (report_diffusion && !my_pe)
         printf("%d var %d sum %lf old %lf diff %lf %lf tol %lf\n",
                ts, v, sum[v], grid_sum[v], (sum[v] - grid_sum[v]),
                (fabs(sum[v] - grid_sum[v])/grid_sum[v]), tol);
      if (stencil || v == 0) {
         if (fabs(sum[v] - grid_sum[v])/grid_sum[v] > tol) {
            if (!my_pe)
               printf("Time step %d sum %lf (old %lf) variable %d difference too large\n", ts, sum[v], grid_sum[v], v);
            MPI_Abort(MPI_COMM_WORLD, -1);
         }
      }
   }
}
