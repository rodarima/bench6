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
// This file is part of miniAMR and is licensed under the terms contained in the ompss-2/LICENSE file.
//
// Copyright (C) 2020-2021 Barcelona Supercomputing Center (BSC)
//

#include <mpi.h>
#include <stdint.h>


/* comm partner variables */
extern double *send_buff, *recv_buff;    /* use in comm and for balancing blocks */
extern double *send_buffers[3];
extern double *recv_buffers[3];

/* comm buffers for exchange */
extern double *send_buff_exchange;
extern double *recv_buff_exchange;
extern int max_exchanges;

extern uint64_t s_buf_size, r_buf_size;

extern int num_comm_partners[3],  // number of comm partners in each direction
    *comm_partner[3],      // list of comm partners in each direction
    max_comm_part[3],      // lengths of comm partners arrays
    *send_size[3],         // send sizes for each comm partner
    *recv_size[3],         // recv sizes for each comm partner
    *comm_index[3],        // index into comm_block, _face_case, and offsets
    *comm_num[3],          // number of blocks for each comm partner
    *comm_block[3],        // array containing local block number for comm
    *comm_face_case[3],    // array containing face cases for comm
    *comm_pos[3],          // position for center of sending face
    *comm_pos1[3],         // perpendicular position of sending face
    num_cases[3],          // amount used in above six arrays
    max_num_cases[3];      // length of above six arrays

extern uint64_t *comm_send_off[3], // offset into send buffer (global, convert to local)
    *comm_recv_off[3],      // offset into recv buffer
    s_buf_num[3],           // total amount being sent in each direction
    r_buf_num[3];           // total amount being received in each direction

extern MPI_Request *request, *s_req;

extern int max_num_req;

// for comm_parent - this is a non-symetric communication

typedef struct {
   int num_comm_part;          // number of other cores to communicate with
   int *comm_part;             // core to communicate with
   int *comm_num;              // number to communicate to each core
   int *index;                 // offset into next two arrays
   num_sz *comm_b;                // block number to communicate from
   num_sz *comm_p;                // parent number of block (for sorting)
   int *comm_c;                // child number of block
   int max_part;               // max communication partners
   int num_cases;              // number to communicate
   int max_cases;              // max number to communicate
} par_comm;
extern par_comm par_b, par_p, par_p1;

extern int *bin, *gbin;

extern MPI_Comm *comms;
extern int *me;
extern int *np;
