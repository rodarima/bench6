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

#include <math.h>
#include <mpi.h>
#include <string.h>

#include "block.h"
#include "comm.h"
#include "proto.h"

void stencil_calc(int, int, int);

// This routine does the stencil calculations.
void stencil_driver(int var, int num, int cacl_stage)
{
   stencil_calc(var, num, stencil);
}

void stencil_calc(int var, int num, int stencil_in)
{
   for (int in = 0; in < sorted_index[num_refine+1]; in++) {
      double *barray = blocks[sorted_list[in].n].array;

      #pragma oss task label("stencil") inout(barray[var*block3d_size;num*block3d_size])
      for (int v = var; v < var+num; v++) {
         typedef double (*block3d_t)[y_block_size+2][z_block_size+2];
         block3d_t matrix = (block3d_t) &(barray[v*block3d_size]);

         double work[x_block_size+2][y_block_size+2][z_block_size+2];
         stencil_kernel((double *) matrix, (double *) work);
      }
   }
   total_fp_divs += (double) num_active*num_cells;
   total_fp_adds += (double) 6*num_active*num_cells;
}
