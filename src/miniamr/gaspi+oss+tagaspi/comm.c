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
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <TAGASPI.h>

#include "block.h"
#include "comm.h"
#include "timer.h"
#include "proto.h"

// The routines in this file are used in the communication of ghost values
// between blocks, both on processor and off processor.

void ack_remote(gaspi_segment_id_t segid, gaspi_notification_id_t notid)
{
    CHECK(tagaspi_notify_async_wait(segid, notid, GASPI_NOTIFICATION_IGNORE));
}

// Main communication routine that sends and recieves ghost values between
// blocks on different processors and directos blocks on the same processor
// to exchange their ghost values.
void comm(int start, int num_comm, int stage, int initial, int final)
{
   int i, j, k, l, m, n, dir, o, in, num, target, rem, g;
   int which, index, type, size, msize;
   int permutations[6][3] = { {0, 1, 2}, {1, 2, 0}, {2, 0, 1},
                              {0, 2, 1}, {1, 0, 2}, {2, 1, 0} };
   uint64_t offset, moffset, roffset;
   double *barray, *barray1;
   double t1, t2, t3, t4;
   block *bp;

   gaspi_notification_id_t lnotid, rnotid;

   const int bstart = start * block3d_size;
   const int bsize = num_comm * block3d_size;

   gaspi_queue_id_t queue = 0;

   for (o = 0; o < 3; o++) {
      if (permute)
         dir = permutations[stage%6][o];
      else
         dir = o;
      type = dir;

      // The default mode is to use a single send/recv buffer for sending/receiving
      // in all directions. In this case, all send_buffers pointers point to send_buff,
      // and all recv_buffers point to recv_buff. In contrast, the separate buffers
      // option makes the application to use different buffers in each direction, so
      // that communication tasks from different directions are concurrent
      double *send_buffer = send_buffers[dir];
      double *recv_buffer = recv_buffers[dir];

      gaspi_segment_id_t send_seg = send_segs[dir];
      gaspi_segment_id_t recv_seg = recv_segs[dir];

      t1 = timer();
      for (i = 0; i < num_comm_partners[dir]; i++) {
         if (send_faces && !max_comm_tasks) { // Single face per message
            for (n = 0; n < comm_num[dir][i]; n++) {
               index = comm_index[dir][i]+n;
               offset = comm_recv_off[dir][index];
               size = comm_face_msg_size(dir, comm_face_case[dir][index], 1);
               lnotid = comm_loc_fid[dir][comm_par_idx[dir][i]+n];

               #pragma oss task label("irecv_face") out(recv_buffer[offset;size])
               {
                  CHECK(tagaspi_notify_async_wait(recv_seg, lnotid, GASPI_NOTIFICATION_IGNORE));
               }
               counter_halo_recv[dir]++;
               size_mesg_recv[dir] += (double) size*sizeof(double);
            }
         } else if (send_faces && max_comm_tasks) { // Multiple faces per message
            target = comm_num[dir][i] / max_comm_tasks;
            rem = comm_num[dir][i] - target * max_comm_tasks;
            uint64_t offsets[target+1];
            uint64_t sizes[target+1];

            n = g = 0;
            while (n < comm_num[dir][i]) {
               num = target + (g < rem);

               moffset = comm_recv_off[dir][comm_index[dir][i]+n];
               lnotid = comm_loc_fid[dir][comm_par_idx[dir][i]+n];
               msize = 0;
               for (m = 0; m < num; m++) {
                  index = comm_index[dir][i]+n+m;
                  offsets[m] = comm_recv_off[dir][index];
                  sizes[m] = comm_face_msg_size(dir, comm_face_case[dir][index], 1);

                  msize += sizes[m];
                  assert(offsets[m] >= moffset && offsets[m] < moffset+msize);

                  size_mesg_recv[dir] += (double) sizes[m]*sizeof(double);
                  counter_halo_recv[dir]++;
               }

               #pragma oss task label("irecv_face") shared(offsets, sizes) \
                     out({recv_buffer[offsets[f];sizes[f]], f=0;num})
               {
                  CHECK(tagaspi_notify_async_wait(recv_seg, lnotid, GASPI_NOTIFICATION_IGNORE));
               }
               n += num; g++;
            }
         } else { // Single message per neighbor
            index = comm_index[dir][i];
            offset = comm_recv_off[dir][index];
            size = recv_size[dir][i];
            lnotid = comm_loc_fid[dir][comm_par_idx[dir][i]];

            #pragma oss task label("irecv_faces") out(recv_buffer[offset;size])
            {
               CHECK(tagaspi_notify_async_wait(recv_seg, lnotid, GASPI_NOTIFICATION_IGNORE));
            }
            counter_halo_recv[dir]++;
            size_mesg_recv[dir] += (double) recv_size[dir][i]*sizeof(double);
         }
      }
      timer_comm_recv[dir] += timer() - t1;

/**** the send and recv list can be same if kept ordered (length can be diff)
        would need to expand face case so that:     **** done ****
        0   whole -> whole
        1   whole -> whole   w/ 27 point stencil (send corners)
        2-5 whole -> quarter w/ number indicating which quarter (matters recv)
        6-9 quarter -> whole w/ number indicating which quarter (matters send)
        + 10 for E, N, U
**** one large send buffer -- can pack and send for a neighbor and reuse */
      for (i = 0; i < num_comm_partners[dir]; i++) {
         if (send_faces && !max_comm_tasks) { // Single message per face
            for (n = 0; n < comm_num[dir][i]; n++) {
               t2 = timer();
               index = comm_index[dir][i]+n;
               offset = comm_send_off[dir][index];
               roffset = comm_rem_off[dir][comm_par_idx[dir][i]+n];
               size = comm_face_msg_size(dir, comm_face_case[dir][index], 0);
               barray = blocks[comm_block[dir][index]].array;
               lnotid = comm_loc_fid[dir][comm_par_idx[dir][i]+n];
               rnotid = comm_rem_fid[dir][comm_par_idx[dir][i]+n];

               #pragma oss task label("pack_face") in(barray[bstart;bsize]) \
                     out(send_buffer[offset;size])
               pack_face(&send_buffer[offset], comm_block[dir][index],
                         comm_face_case[dir][index], dir, start,
                         num_comm);
               counter_face_send[dir]++;
               t3 = timer();

               if (initial) {
                  #pragma oss task label("isend_face") in(send_buffer[offset;size])
                  {
                     CHECK(tagaspi_write_notify(send_seg, offset*sizeof(double),
                           comm_partner[dir][i], recv_seg, roffset*sizeof(double),
                           size*sizeof(double), rnotid, 1, queue));
                  }
               } else {
                  #pragma oss task label("isend_face") in(send_buffer[offset;size]) onready(ack_remote(send_seg, lnotid))
                  {
                     CHECK(tagaspi_write_notify(send_seg, offset*sizeof(double),
                           comm_partner[dir][i], recv_seg, roffset*sizeof(double),
                           size*sizeof(double), rnotid, 1, queue));
                  }
               }
               queue = (queue < num_queues-1) ? queue+1 : 0;

               counter_halo_send[dir]++;
               size_mesg_send[dir] += (double) size*sizeof(double);
               t4 = timer();
               timer_comm_pack[dir] += t3 - t2;
               timer_comm_send[dir] += t4 - t3;
            }
         } else if (send_faces && max_comm_tasks) { // Multiple faces per message
            target = comm_num[dir][i] / max_comm_tasks;
            rem = comm_num[dir][i] - target * max_comm_tasks;
            uint64_t offsets[target+1];
            uint64_t sizes[target+1];

            n = g = 0;
            while (n < comm_num[dir][i]) {
               num = target + (g < rem);

               moffset = comm_send_off[dir][comm_index[dir][i]+n];
               roffset = comm_rem_off[dir][comm_par_idx[dir][i]+n];
               lnotid = comm_loc_fid[dir][comm_par_idx[dir][i]+n];
               rnotid = comm_rem_fid[dir][comm_par_idx[dir][i]+n];
               msize = 0;
               for (m = 0; m < num; m++) {
                  index = comm_index[dir][i]+n+m;
                  offset = offsets[m] = comm_send_off[dir][index];
                  size = sizes[m] = comm_face_msg_size(dir, comm_face_case[dir][index], 0);
                  barray = blocks[comm_block[dir][index]].array;

                  msize += size;
                  assert(offsets[m] >= moffset && offsets[m] < moffset+msize);

                  #pragma oss task label("pack_face") in(barray[bstart;bsize]) \
                        out(send_buffer[offset;size])
                  pack_face(&send_buffer[offset], comm_block[dir][index],
                            comm_face_case[dir][index], dir, start, num_comm);

                  counter_face_send[dir]++;
                  counter_halo_send[dir]++;
                  size_mesg_send[dir] += (double) size*sizeof(double);
               }

               if (initial) {
                  #pragma oss task label("isend_face") shared(offsets, sizes) \
                        in({send_buffer[offsets[f];sizes[f]], f=0;num})
                  {
                     CHECK(tagaspi_write_notify(send_seg, moffset*sizeof(double),
                           comm_partner[dir][i], recv_seg, roffset*sizeof(double),
                           msize*sizeof(double), rnotid, 1, queue));
                  }
               } else {
                  #pragma oss task label("isend_face") shared(offsets, sizes) \
                        in({send_buffer[offsets[f];sizes[f]], f=0;num}) onready(ack_remote(send_seg, lnotid))
                  {
                     CHECK(tagaspi_write_notify(send_seg, moffset*sizeof(double),
                           comm_partner[dir][i], recv_seg, roffset*sizeof(double),
                           msize*sizeof(double), rnotid, 1, queue));
                  }
               }
               queue = (queue < num_queues-1) ? queue+1 : 0;
               n += num; g++;
            }
         } else { // Single message per neighbor
            t2 = timer();
            for (n = 0; n < comm_num[dir][i]; n++) {
               index = comm_index[dir][i]+n;
               offset = comm_send_off[dir][index];
               size = comm_face_msg_size(dir, comm_face_case[dir][index], 0);
               barray = blocks[comm_block[dir][index]].array;

               #pragma oss task label("pack_face") in(barray[bstart;bsize]) \
                     out(send_buffer[offset;size])
               pack_face(&send_buffer[offset], comm_block[dir][index],
                         comm_face_case[dir][index], dir, start, num_comm);
            }
            counter_face_send[dir] += comm_num[dir][i];
            t3 = timer();

            size = send_size[dir][i];
            offset = comm_send_off[dir][comm_index[dir][i]];
            roffset = comm_rem_off[dir][comm_par_idx[dir][i]];
            lnotid = comm_loc_fid[dir][comm_par_idx[dir][i]];
            rnotid = comm_rem_fid[dir][comm_par_idx[dir][i]];
            num = comm_num[dir][i];

            if (initial) {
               #pragma oss task label("isend_faces") in(send_buffer[offset;size])
               {
                  CHECK(tagaspi_write_notify(send_seg, offset*sizeof(double),
                        comm_partner[dir][i], recv_seg, roffset*sizeof(double),
                        size*sizeof(double), rnotid, 1, queue));
               }
            } else {
               #pragma oss task label("isend_faces") in(send_buffer[offset;size]) onready(ack_remote(send_seg, lnotid))
               {
                  CHECK(tagaspi_write_notify(send_seg, offset*sizeof(double),
                        comm_partner[dir][i], recv_seg, roffset*sizeof(double),
                        size*sizeof(double), rnotid, 1, queue));
               }
            }
            queue = (queue < num_queues-1) ? queue+1 : 0;

            counter_halo_send[dir]++;
            size_mesg_send[dir] += (double) size*sizeof(double);
            t4 = timer();
            timer_comm_pack[dir] += t3 - t2;
            timer_comm_send[dir] += t4 - t3;
         }
      }

      // While values are being sent over the mesh, go through and direct
      // blocks to exchange ghost values with other blocks that are on
      // processor.  Also apply boundary conditions for boundary of domain.
      for (in = 0; in < sorted_index[num_refine+1]; in++) {
         bp = &blocks[n = sorted_list[in].n];
         barray = bp->array;
         for (l = dir*2; l < (dir*2 + 2); l++) {
            if (bp->nei_level[l] == bp->level) {
               t2 = timer();
               if ((m = bp->nei[l][0][0]) > n) {
                  barray1 = blocks[m].array;
                  #pragma oss task label("on_proc_comm") \
                        inout(barray[bstart;bsize])      \
                        inout(barray1[bstart;bsize])
                  on_proc_comm(n, m, l, start, num_comm);
                  counter_same[dir] += 2;
               }
               timer_comm_same[dir] += timer() - t2;
            } else if (bp->nei_level[l] == (bp->level+1)) {
               t2 = timer();
               for (i = 0; i < 2; i++)
                  for (j = 0; j < 2; j++)
                     if ((m = bp->nei[l][i][j]) > n) {
                        barray1 = blocks[m].array;
                        #pragma oss task label("on_proc_comm_diff1") \
                              inout(barray[bstart;bsize])            \
                              inout(barray1[bstart;bsize])
                        on_proc_comm_diff(n, m, l, i, j, start, num_comm);
                        counter_diff[dir] += 2;
                     }
               timer_comm_diff[dir] += timer() - t2;
            } else if (bp->nei_level[l] == (bp->level-1)) {
               t2 = timer();
               if ((m = bp->nei[l][0][0]) > n) {
                  k = dir*2 + 1 - l%2;
                  for (i = 0; i < 2; i++)
                     for (j = 0; j < 2; j++)
                        if (blocks[m].nei[k][i][j] == n) {
                           barray1 = blocks[m].array;
                           #pragma oss task label("on_proc_comm_diff2") \
                                 inout(barray[bstart;bsize])            \
                                 inout(barray1[bstart;bsize])
                           on_proc_comm_diff(m, n, k, i, j, start, num_comm);
                           counter_diff[dir] += 2;
                        }
               }
               timer_comm_diff[dir] += timer() - t2;
            } else if (bp->nei_level[l] == -2) {
               t2 = timer();
               #pragma oss task label("boundary_cond") inout(barray[bstart;bsize])
               apply_bc(l, bp, start, num_comm);
               counter_bc[dir]++;
               timer_comm_bc[dir] += timer() - t2;
            } else {
               printf("ERROR: misconnected block\n");
               exit(-1);
            }
         }
      }

      t3 = timer();

      for (i = 0; i < num_comm_partners[dir]; i++) {
         if (send_faces && !max_comm_tasks) {
            for (n = 0; n < comm_num[dir][i]; n++) {
               index = comm_index[dir][i]+n;
               offset = comm_recv_off[dir][index];
               size = comm_face_msg_size(dir, comm_face_case[dir][index], 1);
               bp = &blocks[comm_block[dir][index]];

               if (final) {
                  #pragma oss task label("unpack_face") in(recv_buffer[offset;size]) \
                        inout(bp->array[bstart;bsize])
                  unpack_face(&recv_buffer[offset], comm_block[dir][index],
                              comm_face_case[dir][index], dir, start, num_comm);
               } else {
                  rnotid = comm_rem_fid[dir][comm_par_idx[dir][i]+n];

                  #pragma oss task label("unpack_face") in(recv_buffer[offset;size]) \
                        inout(bp->array[bstart;bsize])
                  {
                     unpack_face(&recv_buffer[offset], comm_block[dir][index],
                                 comm_face_case[dir][index], dir, start, num_comm);
                     CHECK(gaspi_notify(send_seg, comm_partner[dir][i], rnotid, 1, queue, GASPI_BLOCK));
                  }
                  queue = (queue < num_queues-1) ? queue+1 : 0;
               }
            }
         } else {
            target = comm_num[dir][i] / max_comm_tasks;
            rem = comm_num[dir][i] - target * max_comm_tasks;
            uint64_t offsets[target+1];
            uint64_t sizes[target+1];

            n = g = 0;
            while (n < comm_num[dir][i]) {
               num = target + (g < rem);

               for (m = 0; m < num; m++) {
                  index = comm_index[dir][i]+n+m;
                  offset = offsets[m] = comm_recv_off[dir][index];
                  size = sizes[m] = comm_face_msg_size(dir, comm_face_case[dir][index], 1);
                  bp = &blocks[comm_block[dir][index]];

                  #pragma oss task label("unpack_face") in(recv_buffer[offset;size]) \
                        inout(bp->array[bstart;bsize])
                  unpack_face(&recv_buffer[offset], comm_block[dir][index],
                              comm_face_case[dir][index], dir, start, num_comm);
               }

               if (!final) {
                  rnotid = comm_rem_fid[dir][comm_par_idx[dir][i]+n];

                  #pragma oss task label("notify_face") shared(offsets, sizes) \
                        inout({recv_buffer[offsets[f];sizes[f]], f=0;num})
                  {
                     CHECK(gaspi_notify(send_seg, comm_partner[dir][i], rnotid, 1, queue, GASPI_BLOCK));
                  }
                  queue = (queue < num_queues-1) ? queue+1 : 0;
               }
               n += num; ++g;
            }
         }
         counter_face_recv[dir] += comm_num[dir][i];
      }
      t4 = timer();
      timer_comm_unpack[dir] += t4 - t3;

      timer_comm_dir[dir] += timer() - t1;
   }
}

// Pack face to send - note different cases for different directions.
void pack_face(double *send_buf, int block_num, int face_case, int dir,
               int start, int num_comm)
{
   int i, j, k, n, m;
   int is, ie, js, je, ks, ke;
   block *bp;

   typedef double (*block3d_t)[y_block_size+2][z_block_size+2];

   bp = &blocks[block_num];

   if (!code) {

      if (dir == 0) {        /* X - East, West */

         /* X directions (East and West) sent first, so just send
            the real values and no ghosts
         */
         if (face_case >= 10) { /* +X - East */
            i = x_block_size;
            face_case = face_case - 10;
         } else                 /* -X - West */
            i = 1;
         if (face_case < 2) {        /* whole face -> whole face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (j = 1; j <= y_block_size; j++)
                  for (k = 1; k <= z_block_size; k++, n++)
                     send_buf[n] = matrix[i][j][k];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* whole face -> quarter face - case does not matter */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (j = 1; j < y_block_size; j += 2)
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i][j  ][k  ] +
                                   matrix[i][j  ][k+1] +
                                   matrix[i][j+1][k  ] +
                                   matrix[i][j+1][k+1];
            }
         } else {                     /* quarter face -> whole face */
            /* four cases - figure out which quarter of face to send */
            if (face_case%2 == 0) {
               js = 1;
               je = y_block_half;
            } else {
               js = y_block_half + 1;
               je = y_block_size;
            }
            if ((face_case/2)%2 == 1) {
               ks = 1;
               ke = z_block_half;
            } else {
               ks = z_block_half + 1;
               ke = z_block_size;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (j = js; j <= je; j++)
                  for (k = ks; k <= ke; k++, n++)
                     send_buf[n] = matrix[i][j][k]/4.0;
            }
         }

      } else if (dir == 1) { /* Y - North, South */

         /* Y directions (North and South) sent second, so send the real values
         */
         if (face_case >= 10) { /* +Y - North */
            j = y_block_size;
            face_case = face_case - 10;
         } else                 /* -Y - South */
            j = 1;
         if (face_case == 0) {        /* whole face -> whole face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 1; i <= x_block_size; i++)
                  for (k = 1; k <= z_block_size; k++, n++)
                     send_buf[n] = matrix[i][j][k];
            }
         } else if (face_case == 1) {
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 0; i <= x_block_size+1; i++)
                  for (k = 1; k <= z_block_size; k++, n++)
                     send_buf[n] = matrix[i][j][k];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* whole face -> quarter face - case does not matter */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 1; i < x_block_size; i += 2)
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i  ][j][k  ] +
                                   matrix[i  ][j][k+1] +
                                   matrix[i+1][j][k  ] +
                                   matrix[i+1][j][k+1];
            }
         } else {                     /* quarter face -> whole face */
            /* four cases - figure out which quarter of face to send */
            if (face_case%2 == 0) {
               is = 1;
               ie = x_block_half;
            } else {
               is = x_block_half + 1;
               ie = x_block_size;
            }
            if ((face_case/2)%2 == 1) {
               ks = 1;
               ke = z_block_half;
            } else {
               ks = z_block_half + 1;
               ke = z_block_size;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = is; i <= ie; i++)
                  for (k = ks; k <= ke; k++, n++)
                     send_buf[n] = matrix[i][j][k]/4.0;
            }
         }

      } else {               /* Z - Up, Down */

         /* Z directions (Up and Down) sent last
         */
         if (face_case >= 10) { /* +Z - Up */
            k = z_block_size;
            face_case = face_case - 10;
         } else                 /* -Z - Down */
            k = 1;
         if (face_case == 0) {        /* whole face -> whole face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 1; i <= x_block_size; i++)
                  for (j = 1; j <= y_block_size; j++, n++)
                     send_buf[n] = matrix[i][j][k];
            }
         } else if (face_case == 1) {
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 0; i <= x_block_size+1; i++)
                  for (j = 0; j <= y_block_size+1; j++, n++)
                     send_buf[n] = matrix[i][j][k];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* whole face -> quarter face - case does not matter */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 1; i < x_block_size; i += 2)
                  for (j = 1; j < y_block_size; j += 2, n++)
                     send_buf[n] = matrix[i  ][j  ][k] +
                                   matrix[i  ][j+1][k] +
                                   matrix[i+1][j  ][k] +
                                   matrix[i+1][j+1][k];
            }
         } else {                     /* quarter face -> whole face */
            /* four cases - figure out which quarter of face to send */
            if (face_case%2 == 0) {
               is = 1;
               ie = x_block_half;
            } else {
               is = x_block_half + 1;
               ie = x_block_size;
            }
            if ((face_case/2)%2 == 1) {
               js = 1;
               je = y_block_half;
            } else {
               js = y_block_half + 1;
               je = y_block_size;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = is; i <= ie; i++)
                  for (j = js; j <= je; j++, n++)
                     send_buf[n] = matrix[i][j][k]/4.0;
            }
         }
      }

   } else if (code == 1) { /* send all ghosts */

      if (dir == 0) {        /* X - East, West */

         if (face_case >= 10) { /* +X - East */
            i = x_block_size;
            face_case = face_case - 10;
         } else                 /* -X - West */
            i = 1;
         if (face_case < 2) {        /* whole face -> whole face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (j = 0; j <= y_block_size+1; j++)
                  for (k = 0; k <= z_block_size+1; k++, n++)
                     send_buf[n] = matrix[i][j][k];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* whole face -> quarter face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               if (face_case%2 == 0) {
                  j = 0;
                  if ((face_case/2)%2 == 1) {
                     k = 0;
                     send_buf[n++] = matrix[i][j][k];
                  }
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i][j][k  ] +
                                   matrix[i][j][k+1];
                  if ((face_case/2)%2 == 0) {
                     k = z_block_size + 1;
                     send_buf[n++] = matrix[i][j][k];
                  }
               }
               for (j = 1; j < y_block_size; j += 2) {
                  if ((face_case/2)%2 == 1) {
                     k = 0;
                     send_buf[n++] = matrix[i][j  ][k] +
                                     matrix[i][j+1][k];
                  }
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i][j  ][k  ] +
                                   matrix[i][j  ][k+1] +
                                   matrix[i][j+1][k  ] +
                                   matrix[i][j+1][k+1];
                  if ((face_case/2)%2 == 0) {
                     k = z_block_size + 1;
                     send_buf[n++] = matrix[i][j  ][k] +
                                     matrix[i][j+1][k];
                  }
               }
               if (face_case%2 == 1) {
                  j = y_block_size + 1;
                  if ((face_case/2)%2 == 1) {
                     k = 0;
                     send_buf[n++] = matrix[i][j][k];
                  }
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i][j][k  ] +
                                   matrix[i][j][k+1];
                  if ((face_case/2)%2 == 0) {
                     k = z_block_size + 1;
                     send_buf[n++] = matrix[i][j][k];
                  }
               }
            }
         } else {                     /* quarter face -> whole face */
            /* four cases - figure out which quarter of face to send */
            if (face_case%2 == 0) {
               js = 0;
               je = y_block_half + 1;
            } else {
               js = y_block_half;
               je = y_block_size + 1;
            }
            if ((face_case/2)%2 == 1) {
               ks = 0;
               ke = z_block_half + 1;
            } else {
               ks = z_block_half;
               ke = z_block_size + 1;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (j = js; j <= je; j++)
                  for (k = ks; k <= ke; k++, n++)
                     send_buf[n] = matrix[i][j][k]/4.0;
            }
         }

      } else if (dir == 1) { /* Y - North, South */

         if (face_case >= 10) { /* +Y - North */
            j = y_block_size;
            face_case = face_case - 10;
         } else                 /* -Y - South */
            j = 1;
         if (face_case < 2) {        /* whole face -> whole face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 0; i <= x_block_size+1; i++)
                  for (k = 0; k <= z_block_size+1; k++, n++)
                     send_buf[n] = matrix[i][j][k];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* whole face -> quarter face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               if (face_case%2 == 0) {
                  i = 0;
                  if ((face_case/2)%2 == 1) {
                     k = 0;
                     send_buf[n++] = matrix[i][j][k];
                  }
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i][j][k  ] +
                                   matrix[i][j][k+1];
                  if ((face_case/2)%2 == 0) {
                     k = z_block_size + 1;
                     send_buf[n++] = matrix[i][j][k];
                  }
               }
               for (i = 1; i < x_block_size; i += 2) {
                  if ((face_case/2)%2 == 1) {
                     k = 0;
                     send_buf[n++] = matrix[i  ][j][k] +
                                     matrix[i+1][j][k];
                  }
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i  ][j][k  ] +
                                   matrix[i  ][j][k+1] +
                                   matrix[i+1][j][k  ] +
                                   matrix[i+1][j][k+1];
                  if ((face_case/2)%2 == 0) {
                     k = z_block_size + 1;
                     send_buf[n++] = matrix[i  ][j][k] +
                                     matrix[i+1][j][k];
                  }
               }
               if (face_case%2 == 1) {
                  i = x_block_size + 1;
                  if ((face_case/2)%2 == 1) {
                     k = 0;
                     send_buf[n++] = matrix[i][j][k];
                  }
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i][j][k  ] +
                                   matrix[i][j][k+1];
                  if ((face_case/2)%2 == 0) {
                     k = z_block_size + 1;
                     send_buf[n++] = matrix[i][j][k];
                  }
               }
            }
         } else {                     /* quarter face -> whole face */
            /* four cases - figure out which quarter of face to send */
            if (face_case%2 == 0) {
               is = 0;
               ie = x_block_half + 1;
            } else {
               is = x_block_half;
               ie = x_block_size + 1;
            }
            if ((face_case/2)%2 == 1) {
               ks = 0;
               ke = z_block_half + 1;
            } else {
               ks = z_block_half;
               ke = z_block_size + 1;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = is; i <= ie; i++)
                  for (k = ks; k <= ke; k++, n++)
                     send_buf[n] = matrix[i][j][k]/4.0;
            }
         }

      } else {               /* Z - Up, Down */

         /* Z directions (Up and Down) sent last
         */
         if (face_case >= 10) { /* +Z - Up */
            k = z_block_size;
            face_case = face_case - 10;
         } else                 /* -Z - Down */
            k = 1;
         if (face_case < 2) {        /* whole face -> whole face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 0; i <= x_block_size+1; i++)
                  for (j = 0; j <= y_block_size+1; j++, n++)
                     send_buf[n] = matrix[i][j][k];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* whole face -> quarter face - case does not matter */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               if (face_case%2 == 0) {
                  i = 0;
                  if ((face_case/2)%2 == 1) {
                     j = 0;
                     send_buf[n++] = matrix[i][j][k];
                  }
                  for (j = 1; j < y_block_size; j += 2, n++)
                     send_buf[n] = matrix[i][j  ][k] +
                                   matrix[i][j+1][k];
                  if ((face_case/2)%2 == 0) {
                     j = y_block_size + 1;
                     send_buf[n++] = matrix[i][j][k];
                  }
               }
               for (i = 1; i < x_block_size; i += 2) {
                  if ((face_case/2)%2 == 1) {
                     j = 0;
                     send_buf[n++] = matrix[i  ][j][k] +
                                     matrix[i+1][j][k];
                  }
                  for (j = 1; j < y_block_size; j += 2, n++)
                     send_buf[n] = matrix[i  ][j  ][k] +
                                   matrix[i  ][j+1][k] +
                                   matrix[i+1][j  ][k] +
                                   matrix[i+1][j+1][k];
                  if ((face_case/2)%2 == 0) {
                     j = y_block_size + 1;
                     send_buf[n++] = matrix[i  ][j][k] +
                                     matrix[i+1][j][k];
                  }
               }
               if (face_case%2 == 1) {
                  i = x_block_size + 1;
                  if ((face_case/2)%2 == 1) {
                     j = 0;
                     send_buf[n++] = matrix[i][j][k];
                  }
                  for (j = 1; j < y_block_size; j += 2, n++)
                     send_buf[n] = matrix[i][j  ][k] +
                                   matrix[i][j+1][k];
                  if ((face_case/2)%2 == 0) {
                     j = y_block_size + 1;
                     send_buf[n++] = matrix[i][j][k];
                  }
               }
            }
         } else {                     /* quarter face -> whole face */
            /* four cases - figure out which quarter of face to send */
            if (face_case%2 == 0) {
               is = 0;
               ie = x_block_half + 1;
            } else {
               is = x_block_half;
               ie = x_block_size + 1;
            }
            if ((face_case/2)%2 == 1) {
               js = 0;
               je = y_block_half + 1;
            } else {
               js = y_block_half;
               je = y_block_size + 1;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = is; i <= ie; i++)
                  for (j = js; j <= je; j++, n++)
                     send_buf[n] = matrix[i][j][k]/4.0;
            }
         }
      }

   } else { /* code == 2 send all ghosts and do all processing on send side */

      if (dir == 0) {        /* X - East, West */

         if (face_case >= 10) { /* +X - East */
            i = x_block_size;
            face_case = face_case - 10;
         } else                 /* -X - West */
            i = 1;
         if (face_case < 2) {        /* whole face -> whole face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (j = 0; j <= y_block_size+1; j++)
                  for (k = 0; k <= z_block_size+1; k++, n++)
                     send_buf[n] = matrix[i][j][k];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* whole face -> quarter face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               if (face_case%2 == 0) {
                  j = 0;
                  if ((face_case/2)%2 == 1) {
                     k = 0;
                     send_buf[n++] = matrix[i][j][k];
                  }
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i][j][k  ] +
                                   matrix[i][j][k+1];
                  if ((face_case/2)%2 == 0) {
                     k = z_block_size + 1;
                     send_buf[n++] = matrix[i][j][k];
                  }
               }
               for (j = 1; j < y_block_size; j += 2) {
                  if ((face_case/2)%2 == 1) {
                     k = 0;
                     send_buf[n++] = matrix[i][j  ][k] +
                                     matrix[i][j+1][k];
                  }
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i][j  ][k  ] +
                                   matrix[i][j  ][k+1] +
                                   matrix[i][j+1][k  ] +
                                   matrix[i][j+1][k+1];
                  if ((face_case/2)%2 == 0) {
                     k = z_block_size + 1;
                     send_buf[n++] = matrix[i][j  ][k] +
                                     matrix[i][j+1][k];
                  }
               }
               if (face_case%2 == 1) {
                  j = y_block_size + 1;
                  if ((face_case/2)%2 == 1) {
                     k = 0;
                     send_buf[n++] = matrix[i][j][k];
                  }
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i][j][k  ] +
                                   matrix[i][j][k+1];
                  if ((face_case/2)%2 == 0) {
                     k = z_block_size + 1;
                     send_buf[n++] = matrix[i][j][k];
                  }
               }
            }
         } else {                     /* quarter face -> whole face */
            /* four cases - figure out which quarter of face to send */
            if (face_case%2 == 0) {
               js = 1;
               je = y_block_half;
            } else {
               js = y_block_half + 1;
               je = y_block_size;
            }
            if ((face_case/2)%2 == 1) {
               ks = 1;
               ke = z_block_half;
            } else {
               ks = z_block_half + 1;
               ke = z_block_size;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               j = js - 1;
               k = ks - 1;
               send_buf[n++] = matrix[i][j][k]/4.0;
               for (k = ks; k <= ke; k++, n+=2)
                  send_buf[n] = send_buf[n+1] = matrix[i][j][k]/4.0;
               k = ke + 1;
               send_buf[n++] = matrix[i][j][k]/4.0;
               for (j = js; j <= je; j++) {
                  k = ks - 1;
                  send_buf[n++] = matrix[i][j][k]/4.0;
                  for (k = ks; k <= ke; k++, n+=2)
                     send_buf[n] = send_buf[n+1] = matrix[i][j][k]/4.0;
                  k = ke + 1;
                  send_buf[n++] = matrix[i][j][k]/4.0;
                  k = ks - 1;
                  send_buf[n++] = matrix[i][j][k]/4.0;
                  for (k = ks; k <= ke; k++, n+=2)
                     send_buf[n] = send_buf[n+1] = matrix[i][j][k]/4.0;
                  k = ke + 1;
                  send_buf[n++] = matrix[i][j][k]/4.0;
               }
               j = je + 1;
               k = ks - 1;
               send_buf[n++] = matrix[i][j][k]/4.0;
               for (k = ks; k <= ke; k++, n+=2)
                  send_buf[n] = send_buf[n+1] = matrix[i][j][k]/4.0;
               k = ke + 1;
               send_buf[n++] = matrix[i][j][k]/4.0;
            }
         }

      } else if (dir == 1) { /* Y - North, South */

         if (face_case >= 10) { /* +Y - North */
            j = y_block_size;
            face_case = face_case - 10;
         } else                 /* -Y - South */
            j = 1;
         if (face_case < 2) {        /* whole face -> whole face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 0; i <= x_block_size+1; i++)
                  for (k = 0; k <= z_block_size+1; k++, n++)
                     send_buf[n] = matrix[i][j][k];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* whole face -> quarter face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               if (face_case%2 == 0) {
                  i = 0;
                  if ((face_case/2)%2 == 1) {
                     k = 0;
                     send_buf[n++] = matrix[i][j][k];
                  }
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i][j][k  ] +
                                   matrix[i][j][k+1];
                  if ((face_case/2)%2 == 0) {
                     k = z_block_size + 1;
                     send_buf[n++] = matrix[i][j][k];
                  }
               }
               for (i = 1; i < x_block_size; i += 2) {
                  if ((face_case/2)%2 == 1) {
                     k = 0;
                     send_buf[n++] = matrix[i  ][j][k] +
                                     matrix[i+1][j][k];
                  }
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i  ][j][k  ] +
                                   matrix[i  ][j][k+1] +
                                   matrix[i+1][j][k  ] +
                                   matrix[i+1][j][k+1];
                  if ((face_case/2)%2 == 0) {
                     k = z_block_size + 1;
                     send_buf[n++] = matrix[i  ][j][k] +
                                     matrix[i+1][j][k];
                  }
               }
               if (face_case%2 == 1) {
                  i = x_block_size + 1;
                  if ((face_case/2)%2 == 1) {
                     k = 0;
                     send_buf[n++] = matrix[i][j][k];
                  }
                  for (k = 1; k < z_block_size; k += 2, n++)
                     send_buf[n] = matrix[i][j][k  ] +
                                   matrix[i][j][k+1];
                  if ((face_case/2)%2 == 0) {
                     k = z_block_size + 1;
                     send_buf[n++] = matrix[i][j][k];
                  }
               }
            }
         } else {                     /* quarter face -> whole face */
            /* four cases - figure out which quarter of face to send */
            if (face_case%2 == 0) {
               is = 1;
               ie = x_block_half;
            } else {
               is = x_block_half + 1;
               ie = x_block_size;
            }
            if ((face_case/2)%2 == 1) {
               ks = 1;
               ke = z_block_half;
            } else {
               ks = z_block_half + 1;
               ke = z_block_size;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               i = is - 1;
               k = ks - 1;
               send_buf[n++] = matrix[i][j][k]/4.0;
               for (k = ks; k <= ke; k++, n+=2)
                  send_buf[n] = send_buf[n+1] = matrix[i][j][k]/4.0;
               k = ke + 1;
               send_buf[n++] = matrix[i][j][k]/4.0;
               for (i = is; i <= ie; i++) {
                  k = ks - 1;
                  send_buf[n++] = matrix[i][j][k]/4.0;
                  for (k = ks; k <= ke; k++, n+=2)
                     send_buf[n] = send_buf[n+1] = matrix[i][j][k]/4.0;
                  k = ke + 1;
                  send_buf[n++] = matrix[i][j][k]/4.0;
                  k = ks - 1;
                  send_buf[n++] = matrix[i][j][k]/4.0;
                  for (k = ks; k <= ke; k++, n+=2)
                     send_buf[n] = send_buf[n+1] = matrix[i][j][k]/4.0;
                  k = ke + 1;
                  send_buf[n++] = matrix[i][j][k]/4.0;
               }
               i = ie + 1;
               k = ks - 1;
               send_buf[n++] = matrix[i][j][k]/4.0;
               for (k = ks; k <= ke; k++, n+=2)
                  send_buf[n] = send_buf[n+1] = matrix[i][j][k]/4.0;
               k = ke + 1;
               send_buf[n++] = matrix[i][j][k]/4.0;
            }
         }

      } else {               /* Z - Up, Down */

         /* Z directions (Up and Down) sent last
         */
         if (face_case >= 10) { /* +Z - Up */
            k = z_block_size;
            face_case = face_case - 10;
         } else                 /* -Z - Down */
            k = 1;
         if (face_case < 2) {        /* whole face -> whole face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 0; i <= x_block_size+1; i++)
                  for (j = 0; j <= y_block_size+1; j++, n++)
                     send_buf[n] = matrix[i][j][k];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* whole face -> quarter face - case does not matter */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               if (face_case%2 == 0) {
                  i = 0;
                  if ((face_case/2)%2 == 1) {
                     j = 0;
                     send_buf[n++] = matrix[i][j][k];
                  }
                  for (j = 1; j < y_block_size; j += 2, n++)
                     send_buf[n] = matrix[i][j  ][k] +
                                   matrix[i][j+1][k];
                  if ((face_case/2)%2 == 0) {
                     j = y_block_size + 1;
                     send_buf[n++] = matrix[i][j][k];
                  }
               }
               for (i = 1; i < x_block_size; i += 2) {
                  if ((face_case/2)%2 == 1) {
                     j = 0;
                     send_buf[n++] = matrix[i  ][j][k] +
                                     matrix[i+1][j][k];
                  }
                  for (j = 1; j < y_block_size; j += 2, n++)
                     send_buf[n] = matrix[i  ][j  ][k] +
                                   matrix[i  ][j+1][k] +
                                   matrix[i+1][j  ][k] +
                                   matrix[i+1][j+1][k];
                  if ((face_case/2)%2 == 0) {
                     j = y_block_size + 1;
                     send_buf[n++] = matrix[i  ][j][k] +
                                     matrix[i+1][j][k];
                  }
               }
               if (face_case%2 == 1) {
                  i = x_block_size + 1;
                  if ((face_case/2)%2 == 1) {
                     j = 0;
                     send_buf[n++] = matrix[i][j][k];
                  }
                  for (j = 1; j < y_block_size; j += 2, n++)
                     send_buf[n] = matrix[i][j  ][k] +
                                   matrix[i][j+1][k];
                  if ((face_case/2)%2 == 0) {
                     j = y_block_size + 1;
                     send_buf[n++] = matrix[i][j][k];
                  }
               }
            }
         } else {                     /* quarter face -> whole face */
            /* four cases - figure out which quarter of face to send */
            if (face_case%2 == 0) {
               is = 1;
               ie = x_block_half;
            } else {
               is = x_block_half + 1;
               ie = x_block_size;
            }
            if ((face_case/2)%2 == 1) {
               js = 1;
               je = y_block_half;
            } else {
               js = y_block_half + 1;
               je = y_block_size;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               i = is - 1;
               j = js - 1;
               send_buf[n++] = matrix[i][j][k]/4.0;
               for (j = js; j <= je; j++, n+=2)
                  send_buf[n] = send_buf[n+1] = matrix[i][j][k]/4.0;
               j = je + 1;
               send_buf[n++] = matrix[i][j][k]/4.0;
               for (i = is; i <= ie; i++) {
                  j = js - 1;
                  send_buf[n++] = matrix[i][j][k]/4.0;
                  for (j = js; j <= je; j++, n+=2)
                     send_buf[n] = send_buf[n+1] = matrix[i][j][k]/4.0;
                  j = je + 1;
                  send_buf[n++] = matrix[i][j][k]/4.0;
                  j = js - 1;
                  send_buf[n++] = matrix[i][j][k]/4.0;
                  for (j = js; j <= je; j++, n+=2)
                     send_buf[n] = send_buf[n+1] = matrix[i][j][k]/4.0;
                  j = je + 1;
                  send_buf[n++] = matrix[i][j][k]/4.0;
               }
               i = ie + 1;
               j = js - 1;
               send_buf[n++] = matrix[i][j][k]/4.0;
               for (j = js; j <= je; j++, n+=2)
                  send_buf[n] = send_buf[n+1] = matrix[i][j][k]/4.0;
               j = je + 1;
               send_buf[n++] = matrix[i][j][k]/4.0;
            }
         }
      }
   }
}

// Unpack ghost values that have been recieved.
// The sense of the face case is reversed since we are receiving what was sent
void unpack_face(double *recv_buf, int block_num, int face_case, int dir,
                 int start, int num_comm)
{
   int i, j, k, n, m;
   int is, ie, js, je, ks, ke;
   block *bp;

   typedef double (*block3d_t)[y_block_size+2][z_block_size+2];

   bp = &blocks[block_num];

   if (!code) {

      if (dir == 0) {        /* X - East, West */

         /* X directions (East and West)
            just recv the real values and no ghosts
            face_case based on send - so reverse
         */
         if (face_case >= 10) { /* +X - from East */
            i = x_block_size + 1;
            face_case = face_case - 10;
         } else                 /* -X - from West */
            i = 0;
         if (face_case < 2) {        /* whole face -> whole face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (j = 1; j <= y_block_size; j++)
                  for (k = 1; k <= z_block_size; k++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* whole face -> quarter face - one case */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (j = 1; j < y_block_size; j += 2)
                  for (k = 1; k < z_block_size; k += 2, n++)
                     matrix[i][j  ][k  ] =
                     matrix[i][j  ][k+1] =
                     matrix[i][j+1][k  ] =
                     matrix[i][j+1][k+1] = recv_buf[n];
            }
         } else {                     /* quarter face -> whole face */
            /* four cases - figure out which quarter of face to recv */
            if (face_case%2 == 0) {
               js = 1;
               je = y_block_half;
            } else {
               js = y_block_half + 1;
               je = y_block_size;
            }
            if ((face_case/2)%2 == 1) {
               ks = 1;
               ke = z_block_half;
            } else {
               ks = z_block_half + 1;
               ke = z_block_size;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (j = js; j <= je; j++)
                  for (k = ks; k <= ke; k++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         }

      } else if (dir == 1) { /* Y - North, South */

         /* Y directions (North and South) sent second, so recv the real values
         */
         if (face_case >= 10) { /* +Y - from North */
            j = y_block_size + 1;
            face_case = face_case - 10;
         } else                 /* -Y - from South */
            j = 0;
         if (face_case == 0) {        /* whole face -> whole face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 1; i <= x_block_size; i++)
                  for (k = 1; k <= z_block_size; k++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         } else if (face_case == 1) {
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 0; i <= x_block_size+1; i++)
                  for (k = 1; k <= z_block_size; k++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* one case - recv into 4 cells per cell sent */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 1; i < x_block_size; i += 2)
                  for (k = 1; k < z_block_size; k += 2, n++)
                     matrix[i  ][j][k  ] =
                     matrix[i  ][j][k+1] =
                     matrix[i+1][j][k  ] =
                     matrix[i+1][j][k+1] = recv_buf[n];
            }
         } else {                     /* quarter face -> whole face */
            /* whole face -> quarter face - determine case */
            if (face_case%2 == 0) {
               is = 1;
               ie = x_block_half;
            } else {
               is = x_block_half + 1;
               ie = x_block_size;
            }
            if ((face_case/2)%2 == 1) {
               ks = 1;
               ke = z_block_half;
            } else {
               ks = z_block_half + 1;
               ke = z_block_size;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = is; i <= ie; i++)
                  for (k = ks; k <= ke; k++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         }

      } else {               /* Z - Up, Down */

         /* Z directions (Up and Down) sent last
         */
         if (face_case >= 10) { /* +Z - from Up */
            k = z_block_size + 1;
            face_case = face_case - 10;
         } else                 /* -Z - from Down */
            k = 0;
         if (face_case == 0) {        /* whole face -> whole face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 1; i <= x_block_size; i++)
                  for (j = 1; j <= y_block_size; j++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         } else if (face_case == 1) {
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 0; i <= x_block_size+1; i++)
                  for (j = 0; j <= y_block_size+1; j++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* one case - receive into 4 cells */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 1; i < x_block_size; i += 2)
                  for (j = 1; j < y_block_size; j += 2, n++)
                     matrix[i  ][j  ][k] =
                     matrix[i  ][j+1][k] =
                     matrix[i+1][j  ][k] =
                     matrix[i+1][j+1][k] = recv_buf[n];
            }
         } else {                     /* quarter face -> whole face */
            /* whole face -> quarter face - determine case */
            if (face_case%2 == 0) {
               is = 1;
               ie = x_block_half;
            } else {
               is = x_block_half + 1;
               ie = x_block_size;
            }
            if ((face_case/2)%2 == 1) {
               js = 1;
               je = y_block_half;
            } else {
               js = y_block_half + 1;
               je = y_block_size;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = is; i <= ie; i++)
                  for (j = js; j <= je; j++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         }
      }

   } else if (code == 1) {  /* send ghosts */

      if (dir == 0) {        /* X - East, West */

         if (face_case >= 10) { /* +X - from East */
            i = x_block_size + 1;
            face_case = face_case - 10;
         } else                 /* -X - from West */
            i = 0;
         if (face_case < 2) {        /* whole face -> whole */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (j = 0; j <= y_block_size+1; j++)
                  for (k = 0; k <= z_block_size+1; k++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* whole face -> quarter face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               j = 0;
               k = 0;
               matrix[i][j][k] = recv_buf[n++];
               for (k = 1; k < z_block_size; k += 2, n++)
                  matrix[i][j][k  ] =
                  matrix[i][j][k+1] = recv_buf[n];
               k = z_block_size + 1;
               matrix[i][j][k] = recv_buf[n++];
               for (j = 1; j < y_block_size; j += 2) {
                  k = 0;
                  matrix[i][j  ][k] =
                  matrix[i][j+1][k] = recv_buf[n++];
                  for (k = 1; k < z_block_size; k += 2, n++)
                     matrix[i][j  ][k  ] =
                     matrix[i][j  ][k+1] =
                     matrix[i][j+1][k  ] =
                     matrix[i][j+1][k+1] = recv_buf[n];
                  k = z_block_size + 1;
                  matrix[i][j  ][k] =
                  matrix[i][j+1][k] = recv_buf[n++];
               }
               j = y_block_size + 1;
               k = 0;
               matrix[i][j][k] = recv_buf[n++];
               for (k = 1; k < z_block_size; k += 2, n++)
                  matrix[i][j][k  ] =
                  matrix[i][j][k+1] = recv_buf[n];
               k = z_block_size + 1;
               matrix[i][j][k] = recv_buf[n++];
            }
         } else {                     /* quarter face -> whole face */
            /* four cases - figure out which quarter of face to recv */
            if (face_case%2 == 0) {
               js = 0;
               je = y_block_half;
            } else {
               js = y_block_half + 1;
               je = y_block_size + 1;
            }
            if ((face_case/2)%2 == 1) {
               ks = 0;
               ke = z_block_half;
            } else {
               ks = z_block_half + 1;
               ke = z_block_size + 1;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (j = js; j <= je; j++)
                  for (k = ks; k <= ke; k++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         }

      } else if (dir == 1) { /* Y - North, South */

         if (face_case >= 10) { /* +Y - from North */
            j = y_block_size + 1;
            face_case = face_case - 10;
         } else                 /* -Y - from South */
            j = 0;
         if (face_case < 2) {        /* whole face -> whole */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 0; i <= x_block_size+1; i++)
                  for (k = 0; k <= z_block_size+1; k++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* whole face -> quarter face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               i = 0;
               k = 0;
               matrix[i][j][k] = recv_buf[n++];
               for (k = 1; k < z_block_size; k += 2, n++)
                  matrix[i][j][k  ] =
                  matrix[i][j][k+1] = recv_buf[n];
               k = z_block_size + 1;
               matrix[i][j][k] = recv_buf[n++];
               for (i = 1; i < x_block_size; i += 2) {
                  k = 0;
                  matrix[i  ][j][k] =
                  matrix[i+1][j][k] = recv_buf[n++];
                  for (k = 1; k < z_block_size; k += 2, n++)
                     matrix[i  ][j][k  ] =
                     matrix[i  ][j][k+1] =
                     matrix[i+1][j][k  ] =
                     matrix[i+1][j][k+1] = recv_buf[n];
                  k = z_block_size + 1;
                  matrix[i  ][j][k] =
                  matrix[i+1][j][k] = recv_buf[n++];
               }
               i = x_block_size + 1;
               k = 0;
               matrix[i][j][k] = recv_buf[n++];
               for (k = 1; k < z_block_size; k += 2, n++)
                  matrix[i][j][k  ] =
                  matrix[i][j][k+1] = recv_buf[n];
               k = z_block_size + 1;
               matrix[i][j][k] = recv_buf[n++];
            }
         } else {                     /* quarter face -> whole face */
            /* whole face -> quarter face - determine case */
            if (face_case%2 == 0) {
               is = 0;
               ie = x_block_half;
            } else {
               is = x_block_half + 1;
               ie = x_block_size + 1;
            }
            if ((face_case/2)%2 == 1) {
               ks = 0;
               ke = z_block_half;
            } else {
               ks = z_block_half + 1;
               ke = z_block_size + 1;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = is; i <= ie; i++)
                  for (k = ks; k <= ke; k++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         }

      } else {               /* Z - Up, Down */

         if (face_case >= 10) { /* +Z - from Up */
            k = z_block_size + 1;
            face_case = face_case - 10;
         } else                 /* -Z - from Down */
            k = 0;
         if (face_case < 2) {        /* whole face -> whole face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 0; i <= x_block_size+1; i++)
                  for (j = 0; j <= y_block_size+1; j++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         } else if (face_case >= 2 && face_case <= 5) {
            /* whole face -> quarter face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               i = 0;
               j = 0;
               matrix[i][j][k] = recv_buf[n++];
               for (j = 1; j < y_block_size; j += 2, n++)
                  matrix[i][j  ][k] =
                  matrix[i][j+1][k] = recv_buf[n];
               j = y_block_size + 1;
               matrix[i][j][k] = recv_buf[n++];
               for (i = 1; i < x_block_size; i += 2) {
                  j = 0;
                  matrix[i  ][j][k] =
                  matrix[i+1][j][k] = recv_buf[n++];
                  for (j = 1; j < y_block_size; j += 2, n++)
                     matrix[i  ][j  ][k] =
                     matrix[i  ][j+1][k] =
                     matrix[i+1][j  ][k] =
                     matrix[i+1][j+1][k] = recv_buf[n];
                  j = y_block_size + 1;
                  matrix[i  ][j][k] =
                  matrix[i+1][j][k] = recv_buf[n++];
               }
               i = x_block_size + 1;
               j = 0;
               matrix[i][j][k] = recv_buf[n++];
               for (j = 1; j < y_block_size; j += 2, n++)
                  matrix[i][j  ][k] =
                  matrix[i][j+1][k] = recv_buf[n];
               j = y_block_size + 1;
               matrix[i][j][k] = recv_buf[n++];
            }
         } else {                     /* quarter face -> whole face */
            /* whole face -> quarter face - determine case */
            if (face_case%2 == 0) {
               is = 0;
               ie = x_block_half;
            } else {
               is = x_block_half + 1;
               ie = x_block_size + 1;
            }
            if ((face_case/2)%2 == 1) {
               js = 0;
               je = y_block_half;
            } else {
               js = y_block_half + 1;
               je = y_block_size + 1;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = is; i <= ie; i++)
                  for (j = js; j <= je; j++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         }
      }

   } else {  /* code == 2 send ghosts and process on send */

      if (dir == 0) {        /* X - East, West */

         if (face_case >= 10) { /* +X - from East */
            i = x_block_size + 1;
            face_case = face_case - 10;
         } else                 /* -X - from West */
            i = 0;
         if (face_case <= 5) {        /* whole face -> whole or quarter face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (j = 0; j <= y_block_size+1; j++)
                  for (k = 0; k <= z_block_size+1; k++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         } else {                     /* quarter face -> whole face */
            /* four cases - figure out which quarter of face to recv */
            if (face_case%2 == 0) {
               js = 0;
               je = y_block_half;
            } else {
               js = y_block_half + 1;
               je = y_block_size + 1;
            }
            if ((face_case/2)%2 == 1) {
               ks = 0;
               ke = z_block_half;
            } else {
               ks = z_block_half + 1;
               ke = z_block_size + 1;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (j = js; j <= je; j++)
                  for (k = ks; k <= ke; k++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         }

      } else if (dir == 1) { /* Y - North, South */

         if (face_case >= 10) { /* +Y - from North */
            j = y_block_size + 1;
            face_case = face_case - 10;
         } else                 /* -Y - from South */
            j = 0;
         if (face_case <= 5) {        /* whole face -> whole or quarter face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 0; i <= x_block_size+1; i++)
                  for (k = 0; k <= z_block_size+1; k++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         } else {                     /* quarter face -> whole face */
            /* whole face -> quarter face - determine case */
            if (face_case%2 == 0) {
               is = 0;
               ie = x_block_half;
            } else {
               is = x_block_half + 1;
               ie = x_block_size + 1;
            }
            if ((face_case/2)%2 == 1) {
               ks = 0;
               ke = z_block_half;
            } else {
               ks = z_block_half + 1;
               ke = z_block_size + 1;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = is; i <= ie; i++)
                  for (k = ks; k <= ke; k++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         }

      } else {               /* Z - Up, Down */

         if (face_case >= 10) { /* +Z - from Up */
            k = z_block_size + 1;
            face_case = face_case - 10;
         } else                 /* -Z - from Down */
            k = 0;
         if (face_case <= 5) {        /* whole face -> whole or quarter face */
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = 0; i <= x_block_size+1; i++)
                  for (j = 0; j <= y_block_size+1; j++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         } else {                     /* quarter face -> whole face */
            /* whole face -> quarter face - determine case */
            if (face_case%2 == 0) {
               is = 0;
               ie = x_block_half;
            } else {
               is = x_block_half + 1;
               ie = x_block_size + 1;
            }
            if ((face_case/2)%2 == 1) {
               js = 0;
               je = y_block_half;
            } else {
               js = y_block_half + 1;
               je = y_block_size + 1;
            }
            for (n = 0, m = start; m < start+num_comm; m++) {
               block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
               for (i = is; i <= ie; i++)
                  for (j = js; j <= je; j++, n++)
                     matrix[i][j][k] = recv_buf[n];
            }
         }
      }
   }
}

// Routine that does on processor communication between two blocks that
// are at the same level of refinement.
void on_proc_comm(int n, int n1, int l, int start, int num_comm)
{
   int i, j, k, m;
   int is, ie, js, je;
   block *bp, *bp1;

   typedef double (*block3d_t)[y_block_size+2][z_block_size+2];

   /* Determine direction and then exchange data across the face
   */
   if (!code) {
      if ((l/2) == 0) {         /* West, East */
         if ((l%2) == 0) {      /* West */
            bp = &blocks[n];
            bp1 = &blocks[n1];
         } else {               /* East */
            bp1 = &blocks[n];
            bp = &blocks[n1];
         }
         copy_ghosts_westeast(bp, bp1, start, num_comm);
      } else if ((l/2) == 1) {  /* South, North */
         if ((l%2) == 0) {      /* South */
            bp = &blocks[n];
            bp1 = &blocks[n1];
         } else {               /* North */
            bp1 = &blocks[n];
            bp = &blocks[n1];
         }
         if (stencil == 7) {
           is = 1;
           ie = x_block_size;
         } else {
           is = 0;
           ie = x_block_size + 1;
         }
         copy_ghosts_southnorth(bp, bp1, start, num_comm, is, ie);
      } else if ((l/2) == 2) {  /* Down, Up */
         if ((l%2) == 0) {      /* Down */
            bp = &blocks[n];
            bp1 = &blocks[n1];
         } else {               /* Up */
            bp1 = &blocks[n];
            bp = &blocks[n1];
         }
         if (stencil == 7) {
           is = 1;
           ie = x_block_size;
           js = 1;
           je = y_block_size;
         } else {
           is = 0;
           ie = x_block_size + 1;
           js = 0;
           je = y_block_size + 1;
         }
         copy_ghosts_downup(bp, bp1, start, num_comm, is, ie, js, je);
      }
   } else {  /* set all ghosts */
      if ((l/2) == 0) {         /* West, East */
         if ((l%2) == 0) {      /* West */
            bp = &blocks[n];
            bp1 = &blocks[n1];
         } else {               /* East */
            bp1 = &blocks[n];
            bp = &blocks[n1];
         }
         for (m = start; m < start+num_comm; m++) {
            block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
            block3d_t matrix1 = (block3d_t) &(bp1->array[m*block3d_size]);
            for (j = 0; j <= y_block_size+1; j++)
               for (k = 0; k <= z_block_size+1; k++) {
                  matrix1[x_block_size+1][j][k] = matrix[1][j][k];
                  matrix[0][j][k] = matrix1[x_block_size][j][k];
               }
         }
      } else if ((l/2) == 1) {  /* South, North */
         if ((l%2) == 0) {      /* South */
            bp = &blocks[n];
            bp1 = &blocks[n1];
         } else {               /* North */
            bp1 = &blocks[n];
            bp = &blocks[n1];
         }
         for (m = start; m < start+num_comm; m++) {
            block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
            block3d_t matrix1 = (block3d_t) &(bp1->array[m*block3d_size]);
            for (i = 0; i <= x_block_size+1; i++)
               for (k = 0; k <= z_block_size+1; k++) {
                  matrix1[i][y_block_size+1][k] = matrix[i][1][k];
                  matrix[i][0][k] = matrix1[i][y_block_size][k];
               }
         }
      } else if ((l/2) == 2) {  /* Down, Up */
         if ((l%2) == 0) {      /* Down */
            bp = &blocks[n];
            bp1 = &blocks[n1];
         } else {               /* Up */
            bp1 = &blocks[n];
            bp = &blocks[n1];
         }
         for (m = start; m < start+num_comm; m++) {
            block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
            block3d_t matrix1 = (block3d_t) &(bp1->array[m*block3d_size]);
            for (i = 0; i <= x_block_size+1; i++)
               for (j = 0; j <= y_block_size+1; j++) {
                  matrix1[i][j][z_block_size+1] = matrix[i][j][1];
                  matrix[i][j][0] = matrix1[i][j][z_block_size];
               }
         }
      }
   }
}

// Routine that does on processor communication between two blocks that are
// at different levels of refinement.  The order of the blocks that are
// being input determine which block is at a higher level of refinement.
void on_proc_comm_diff(int n, int n1, int l, int iq, int jq,
                       int start, int num_comm)
{
   int i, j, k, m;
   int i0, i1, i2, i3, j0, j1, j2, j3, k0, k1, k2, k3;
   block *bp, *bp1;

   typedef double (*block3d_t)[y_block_size+2][z_block_size+2];

   bp = &blocks[n];
   bp1 = &blocks[n1];

   /* (iq, jq) quarter face on block n to whole face on block n1
   */
   if (!code) {
      /* only have to communicate ghost values - bp is level, bp1 is level+1 -
       * in 2 to 1 case get 0..block_half from one proc and
       *                block_half+1..block_size+1 from another
       * in 1 to 2 case get 0..block_size+1 from 0..block_half+1 or
       *                block_half..block_size+1 with averages
       */
      if ((l/2) == 0) {
         if (l == 0) {             /* West */
            i0 = 0;
            i1 = 1;
            i2 = x_block_size + 1;
            i3 = x_block_size;
         } else {                  /* East */
            i0 = x_block_size + 1;
            i1 = x_block_size;
            i2 = 0;
            i3 = 1;
         }
         j1 = jq*y_block_half;
         k1 = iq*z_block_half;
         for (m = start; m < start+num_comm; m++) {
            block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
            block3d_t matrix1 = (block3d_t) &(bp1->array[m*block3d_size]);
            for (j = 1; j <= y_block_half; j++)
               for (k = 1; k <= z_block_half; k++) {
                  matrix1[i2][2*j-1][2*k-1] =
                  matrix1[i2][2*j-1][2*k  ] =
                  matrix1[i2][2*j  ][2*k-1] =
                  matrix1[i2][2*j  ][2*k  ] =
                                             matrix[i1][j+j1][k+k1]/4.0;
                  matrix[i0][j+j1][k+k1] =
                                             matrix1[i3][2*j-1][2*k-1] +
                                             matrix1[i3][2*j-1][2*k  ] +
                                             matrix1[i3][2*j  ][2*k-1] +
                                             matrix1[i3][2*j  ][2*k  ];
               }
         }
      } else if ((l/2) == 1) {
         if (l == 2) {             /* South */
            j0 = 0;
            j1 = 1;
            j2 = y_block_size + 1;
            j3 = y_block_size;
         } else {                  /* North */
            j0 = y_block_size + 1;
            j1 = y_block_size;
            j2 = 0;
            j3 = 1;
         }
         i1 = jq*x_block_half;
         k1 = iq*z_block_half;
         for (m = start; m < start+num_comm; m++) {
            block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
            block3d_t matrix1 = (block3d_t) &(bp1->array[m*block3d_size]);
            for (i = 1; i <= x_block_half; i++)
               for (k = 1; k <= z_block_half; k++) {
                  matrix1[2*i-1][j2][2*k-1] =
                  matrix1[2*i-1][j2][2*k  ] =
                  matrix1[2*i  ][j2][2*k-1] =
                  matrix1[2*i  ][j2][2*k  ] =
                                             matrix[i+i1][j1][k+k1]/4.0;
                  matrix[i+i1][j0][k+k1] =
                                             matrix1[2*i-1][j3][2*k-1] +
                                             matrix1[2*i-1][j3][2*k  ] +
                                             matrix1[2*i  ][j3][2*k-1] +
                                             matrix1[2*i  ][j3][2*k  ];
               }
         }
      } else if ((l/2) == 2) {
         if (l == 4) {             /* Down */
            k0 = 0;
            k1 = 1;
            k2 = z_block_size + 1;
            k3 = z_block_size;
         } else {                  /* Up */
            k0 = z_block_size + 1;
            k1 = z_block_size;
            k2 = 0;
            k3 = 1;
         }
         i1 = jq*x_block_half;
         j1 = iq*y_block_half;
         for (m = start; m < start+num_comm; m++) {
            block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
            block3d_t matrix1 = (block3d_t) &(bp1->array[m*block3d_size]);
            for (i = 1; i <= x_block_half; i++)
               for (j = 1; j <= y_block_half; j++) {
                  matrix1[2*i-1][2*j-1][k2] =
                  matrix1[2*i-1][2*j  ][k2] =
                  matrix1[2*i  ][2*j-1][k2] =
                  matrix1[2*i  ][2*j  ][k2] =
                                              matrix[i+i1][j+j1][k1]/4.0;
                  matrix[i+i1][j+j1][k0] =
                                              matrix1[2*i-1][2*j-1][k3] +
                                              matrix1[2*i-1][2*j  ][k3] +
                                              matrix1[2*i  ][2*j-1][k3] +
                                              matrix1[2*i  ][2*j  ][k3];
               }
         }
      }
   } else {  /* transfer ghosts */
      if ((l/2) == 0) {
         if (l == 0) {             /* West */
            i0 = 0;
            i1 = 1;
            i2 = x_block_size + 1;
            i3 = x_block_size;
         } else {                  /* East */
            i0 = x_block_size + 1;
            i1 = x_block_size;
            i2 = 0;
            i3 = 1;
         }
         j1 = jq*y_block_half;
         k1 = iq*z_block_half;
         j2 = y_block_size + 1;
         j3 = y_block_half + 1;
         k2 = z_block_size + 1;
         k3 = z_block_half + 1;
         for (m = start; m < start+num_comm; m++) {
            block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
            block3d_t matrix1 = (block3d_t) &(bp1->array[m*block3d_size]);

            matrix1[i2][0][0] = matrix[i1][j1][k1]/4.0;
            for (k = 1; k <= z_block_half; k++)
               matrix1[i2][0][2*k-1] =
               matrix1[i2][0][2*k  ] = matrix[i1][j1][k+k1]/4.0;
            matrix1[i2][0][k2] = matrix[i1][j1][k3+k1]/4.0;
            if (jq == 0) {
               if (iq == 0)
                  matrix[i0][0][0 ] = matrix1[i3][0][0 ];
               else
                  matrix[i0][0][k2] = matrix1[i3][0][k2];
               for (k = 1; k <= z_block_half; k++)
                  matrix[i0][0][k+k1] = (matrix1[i3][0][2*k-1] +
                                               matrix1[i3][0][2*k  ]);
            }
            for (j = 1; j <= y_block_half; j++) {
               matrix1[i2][2*j-1][0] =
               matrix1[i2][2*j  ][0] = matrix[i1][j+j1][k1]/4.0;
               if (iq == 0)
                  matrix[i0][j+j1][0 ] = (matrix1[i3][2*j-1][0 ] +
                                                matrix1[i3][2*j  ][0 ]);
               else
                  matrix[i0][j+j1][k2] = (matrix1[i3][2*j-1][k2] +
                                                matrix1[i3][2*j  ][k2]);
               for (k = 1; k <= z_block_half; k++) {
                  matrix1[i2][2*j-1][2*k-1] =
                  matrix1[i2][2*j-1][2*k  ] =
                  matrix1[i2][2*j  ][2*k-1] =
                  matrix1[i2][2*j  ][2*k  ] =
                                             matrix[i1][j+j1][k+k1]/4.0;
                  matrix[i0][j+j1][k+k1] =
                                             matrix1[i3][2*j-1][2*k-1] +
                                             matrix1[i3][2*j-1][2*k  ] +
                                             matrix1[i3][2*j  ][2*k-1] +
                                             matrix1[i3][2*j  ][2*k  ];
               }
               matrix1[i2][2*j-1][k2] =
               matrix1[i2][2*j  ][k2] = matrix[i1][j+j1][k3+k1]/4.0;
            }
            matrix1[i2][j2][0] = matrix[i1][j3+j1][k1]/4.0;
            for (k = 1; k <= z_block_half; k++)
               matrix1[i2][j2][2*k-1] =
               matrix1[i2][j2][2*k  ] = matrix[i1][j3+j1][k+k1]/4.0;
            matrix1[i2][j2][k2] = matrix[i1][j3+j1][k3+k1]/4.0;
            if (jq == 1) {
               if (iq == 0)
                  matrix[i0][j2][0 ] = matrix1[i3][j2][0 ];
               else
                  matrix[i0][j2][k2] = matrix1[i3][j2][k2];
               for (k = 1; k <= z_block_half; k++)
                  matrix[i0][j2][k+k1] = (matrix1[i3][j2][2*k-1] +
                                                matrix1[i3][j2][2*k  ]);
            }
         }
      } else if ((l/2) == 1) {
         if (l == 2) {             /* South */
            j0 = 0;
            j1 = 1;
            j2 = y_block_size + 1;
            j3 = y_block_size;
         } else {                  /* North */
            j0 = y_block_size + 1;
            j1 = y_block_size;
            j2 = 0;
            j3 = 1;
         }
         i1 = jq*x_block_half;
         k1 = iq*z_block_half;
         i2 = x_block_size + 1;
         i3 = x_block_half + 1;
         k2 = z_block_size + 1;
         k3 = z_block_half + 1;
         for (m = start; m < start+num_comm; m++) {
            block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
            block3d_t matrix1 = (block3d_t) &(bp1->array[m*block3d_size]);

            matrix1[0][j2][0 ] = matrix[i1][j1][k1]/4.0;
            for (k = 1; k <= z_block_half; k++)
               matrix1[0][j2][2*k-1] =
               matrix1[0][j2][2*k  ] = matrix[i1][j1][k+k1]/4.0;
            matrix1[0][j2][k2] = matrix[i1][j1][k3+k1]/4.0;
            if (jq == 0) {
               if (iq == 0)
                  matrix[0][j0][0 ] = matrix1[0][j3][0 ];
               else
                  matrix[0][j0][k2] = matrix1[0][j3][k2];
               for (k = 1; k <= z_block_half; k++)
                  matrix[0][j0][k+k1] = (matrix1[0][j3][2*k-1] +
                                               matrix1[0][j3][2*k  ]);
            }
            for (i = 1; i <= x_block_half; i++) {
               matrix1[2*i-1][j2][0] =
               matrix1[2*i  ][j2][0] = matrix[i+i1][j1][k1]/4.0;
               if (iq == 0)
                  matrix[i+i1][j0][0 ] = (matrix1[2*i-1][j3][0 ] +
                                                matrix1[2*i  ][j3][0 ]);
               else
                  matrix[i+i1][j0][k2] = (matrix1[2*i-1][j3][k2] +
                                                matrix1[2*i  ][j3][k2]);
               for (k = 1; k <= z_block_half; k++) {
                  matrix1[2*i-1][j2][2*k-1] =
                  matrix1[2*i-1][j2][2*k  ] =
                  matrix1[2*i  ][j2][2*k-1] =
                  matrix1[2*i  ][j2][2*k  ] =
                                             matrix[i+i1][j1][k+k1]/4.0;
                  matrix[i+i1][j0][k+k1] =
                                             matrix1[2*i-1][j3][2*k-1] +
                                             matrix1[2*i-1][j3][2*k  ] +
                                             matrix1[2*i  ][j3][2*k-1] +
                                             matrix1[2*i  ][j3][2*k  ];
               }
               matrix1[2*i-1][j2][k2] =
               matrix1[2*i  ][j2][k2] = matrix[i+i1][j1][k3+k1]/4.0;
            }
            matrix1[i2][j2][0 ] = matrix[i3+i1][j1][k1]/4.0;
            for (k = 1; k <= z_block_half; k++)
               matrix1[i2][j2][2*k-1] =
               matrix1[i2][j2][2*k  ] = matrix[i3+i1][j1][k+k1]/4.0;
            matrix1[i2][j2][k2] = matrix[i3+i1][j1][k3+k1]/4.0;
            if (jq == 1) {
               if (iq == 0)
                  matrix[i2][j0][0 ] = matrix1[i2][j3][0 ];
               else
                  matrix[i2][j0][k2] = matrix1[i2][j3][k2];
               for (k = 1; k <= z_block_half; k++)
                  matrix[i2][j0][k+k1] = (matrix1[i2][j3][2*k-1] +
                                                matrix1[i2][j3][2*k  ]);
            }
         }
      } else if ((l/2) == 2) {
         if (l == 4) {             /* Down */
            k0 = 0;
            k1 = 1;
            k2 = z_block_size + 1;
            k3 = z_block_size;
         } else {                  /* Up */
            k0 = z_block_size + 1;
            k1 = z_block_size;
            k2 = 0;
            k3 = 1;
         }
         i1 = jq*x_block_half;
         j1 = iq*y_block_half;
         i2 = x_block_size + 1;
         i3 = x_block_half + 1;
         j2 = y_block_size + 1;
         j3 = y_block_half + 1;
         for (m = start; m < start+num_comm; m++) {
            block3d_t matrix = (block3d_t) &(bp->array[m*block3d_size]);
            block3d_t matrix1 = (block3d_t) &(bp1->array[m*block3d_size]);

            matrix1[0][0 ][k2] = matrix[i1][j1][k1]/4.0;
            for (j = 1; j <= y_block_half; j++)
               matrix1[0][2*j-1][k2] =
               matrix1[0][2*j  ][k2] = matrix[i1][j+j1][k1]/4.0;
            matrix1[0][j2][k2] = matrix[i1][j3+j1][k1]/4.0;
            if (jq == 0) {
               if (iq == 0)
                  matrix[0][0 ][k0] = matrix1[0][0 ][k3];
               else
                  matrix[0][j2][k0] = matrix1[0][j2][k3];
               for (j = 1; j <= y_block_half; j++)
                  matrix[0][j+j1][k0] = (matrix1[0][2*j-1][k3] +
                                               matrix1[0][2*j  ][k3]);
            }
            for (i = 1; i <= x_block_half; i++) {
               matrix1[2*i-1][0][k2] =
               matrix1[2*i  ][0][k2] = matrix[i+i1][j1][k1]/4.0;
               if (iq == 0)
                  matrix[i+i1][0][k0] = (matrix1[2*i-1][0][k3] +
                                               matrix1[2*i  ][0][k3]);
               else
                  matrix[i+i1][j2][k0] = (matrix1[2*i-1][j2][k3] +
                                                matrix1[2*i  ][j2][k3]);
               for (j = 1; j <= y_block_half; j++) {
                  matrix1[2*i-1][2*j-1][k2] =
                  matrix1[2*i-1][2*j  ][k2] =
                  matrix1[2*i  ][2*j-1][k2] =
                  matrix1[2*i  ][2*j  ][k2] =
                                              matrix[i+i1][j+j1][k1]/4.0;
                  matrix[i+i1][j+j1][k0] =
                                              matrix1[2*i-1][2*j-1][k3] +
                                              matrix1[2*i-1][2*j  ][k3] +
                                              matrix1[2*i  ][2*j-1][k3] +
                                              matrix1[2*i  ][2*j  ][k3];
               }
               matrix1[2*i-1][j2][k2] =
               matrix1[2*i  ][j2][k2] = matrix[i+i1][j3+j1][k1]/4.0;
            }
            matrix1[i2][0 ][k2] = matrix[i3+i1][j1][k1]/4.0;
            for (j = 1; j <= y_block_half; j++)
               matrix1[i2][2*j-1][k2] =
               matrix1[i2][2*j  ][k2] = matrix[i3+i1][j+j1][k1]/4.0;
            matrix1[i2][j2][k2] = matrix[i3+i1][j3+j1][k1]/4.0;
            if (jq == 1) {
               if (iq == 0)
                  matrix[i2][0 ][k0] = matrix1[i2][0 ][k3];
               else
                  matrix[i2][j2][k0] = matrix1[i2][j2][k3];
               for (j = 1; j <= y_block_half; j++)
                  matrix[i2][j+j1][k0] = (matrix1[i2][2*j-1][k3] +
                                                matrix1[i2][2*j  ][k3]);
            }
         }
      }
   }
}

// Apply reflective boundary conditions to a face of a block.
void apply_bc(int l, block *bp, int start, int num_comm)
{
   int var, i, j, k, f, t;

   typedef double (*block3d_t)[y_block_size+2][z_block_size+2];

   t = 0;
   f = 1;
   if (!code && stencil == 7)
      switch (l) {
         case 1: t = x_block_size + 1;
                 f = x_block_size;
         case 0: for (var = start; var < start+num_comm; var++) {
                    block3d_t matrix = (block3d_t) &(bp->array[var*block3d_size]);
                    for (j = 1; j <= y_block_size; j++)
                       for (k = 1; k <= z_block_size; k++)
                          matrix[t][j][k] = matrix[f][j][k];
                 }
                 break;
         case 3: t = y_block_size + 1;
                 f = y_block_size;
         case 2: for (var = start; var < start+num_comm; var++) {
                    block3d_t matrix = (block3d_t) &(bp->array[var*block3d_size]);
                    for (i = 1; i <= x_block_size; i++)
                       for (k = 1; k <= z_block_size; k++)
                          matrix[i][t][k] = matrix[i][f][k];
                 }
                 break;
         case 5: t = z_block_size + 1;
                 f = z_block_size;
         case 4: for (var = start; var < start+num_comm; var++) {
                    block3d_t matrix = (block3d_t) &(bp->array[var*block3d_size]);
                    for (i = 1; i <= x_block_size; i++)
                       for (j = 1; j <= y_block_size; j++)
                          matrix[i][j][t] = matrix[i][j][f];
                 }
                 break;
      }
   else
      switch (l) {
         case 1: t = x_block_size + 1;
                 f = x_block_size;
         case 0: for (var = start; var < start+num_comm; var++) {
                    block3d_t matrix = (block3d_t) &(bp->array[var*block3d_size]);
                    for (j = 0; j <= y_block_size+1; j++)
                       for (k = 0; k <= z_block_size+1; k++)
                          matrix[t][j][k] = matrix[f][j][k];
                 }
                 break;
         case 3: t = y_block_size + 1;
                 f = y_block_size;
         case 2: for (var = start; var < start+num_comm; var++) {
                    block3d_t matrix = (block3d_t) &(bp->array[var*block3d_size]);
                    for (i = 0; i <= x_block_size+1; i++)
                       for (k = 0; k <= z_block_size+1; k++)
                          matrix[i][t][k] = matrix[i][f][k];
                 }
                 break;
         case 5: t = z_block_size + 1;
                 f = z_block_size;
         case 4: for (var = start; var < start+num_comm; var++) {
                    block3d_t matrix = (block3d_t) &(bp->array[var*block3d_size]);
                    for (i = 0; i <= x_block_size+1; i++)
                       for (j = 0; j <= y_block_size+1; j++)
                          matrix[i][j][t] = matrix[i][j][f];
                 }
                 break;
      }
}
