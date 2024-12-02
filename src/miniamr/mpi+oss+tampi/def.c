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

#include "block.h"
#include "comm.h"
#include "timer.h"

// block.h
block *blocks;
parent *parents;
sorted_block *sorted_list;
int *sorted_index;
int my_pe;
int num_pes;
int max_num_blocks;
int num_refine;
int uniform_refine;
int x_block_size, y_block_size, z_block_size;
int block3d_size;
int num_cells;
int num_vars;
int mat;
int comm_vars;
int init_block_x, init_block_y, init_block_z;
int reorder;
int npx, npy, npz;
int inbalance;
int refine_freq;
int report_diffusion;
int checksum_freq;
int checksum_opt;
int stages_per_ts;
int error_tol;
int num_tsteps;
int stencil;
int report_perf;
int plot_freq;
int lb_opt;
int lb_freq;
int block_change;
int code;
int permute;
int nonblocking;
int refine_ghost;
int use_time;
double end_time;
int send_faces;
int max_comm_tasks;
int separate_buffers;
int change_dir;
int group_blocks;
int limit_move;
int first;
int *dirs;
int max_num_parents;
int num_parents;
int max_active_parent;
int cur_max_level;
num_sz *num_blocks;
num_sz *local_num_blocks;
num_sz *block_start;
int num_active;
int max_active_block;
num_sz global_active;
int x_block_half, y_block_half, z_block_half;
double tol;
double *grid_sum;
int *p8, *p2;
int mesh_size[3];
int max_mesh_size;
int *from, *to;
int msg_len[3][4];
int local_max_b;
int global_max_b;
int lb_countdown;
double *alpha, beta;
double total_fp_divs;
double total_fp_adds;
double total_fp_muls;
int num_objects;
object *objects;
int num_dots;
int max_num_dots;
int max_active_dot;
dot *dots;

// comm.h
double *send_buff, *recv_buff;
double *send_buffers[3];
double *recv_buffers[3];
double *send_buff_exchange;
double *recv_buff_exchange;
int max_exchanges;
uint64_t s_buf_size, r_buf_size;
int num_comm_partners[3],
    *comm_partner[3],
    max_comm_part[3],
    *send_size[3],
    *recv_size[3],
    *comm_index[3],
    *comm_num[3],
    *comm_block[3],
    *comm_face_case[3],
    *comm_pos[3],
    *comm_pos1[3],
    num_cases[3],
    max_num_cases[3];
uint64_t *comm_send_off[3],
    *comm_recv_off[3],
    s_buf_num[3],
    r_buf_num[3];
MPI_Request *request, *s_req;
int max_num_req;
par_comm par_b, par_p, par_p1;
int *bin, *gbin;
MPI_Comm *comms;
int *me;
int *np;

// timer.h
double average[146];
double stddev[143];
double minimum[143];
double maximum[143];
double timer_all;
double timer_comm_all;
double timer_comm_dir[3];
double timer_comm_recv[3];
double timer_comm_pack[3];
double timer_comm_send[3];
double timer_comm_same[3];
double timer_comm_diff[3];
double timer_comm_bc[3];
double timer_comm_wait[3];
double timer_comm_unpack[3];
double timer_calc_all;
double timer_cs_all;
double timer_cs_red;
double timer_cs_calc;
double timer_refine_all;
double timer_refine_co;
double timer_refine_mr;
double timer_refine_cc;
double timer_refine_sb;
double timer_refine_c1;
double timer_refine_c2;
double timer_refine_sy;
double timer_cb_all;
double timer_cb_cb;
double timer_cb_pa;
double timer_cb_mv;
double timer_cb_un;
double timer_lb_all;
double timer_lb_sort;
double timer_lb_pa;
double timer_lb_mv;
double timer_lb_un;
double timer_lb_misc;
double timer_lb_mb;
double timer_lb_ma;
double timer_rs_all;
double timer_rs_ca;
double timer_rs_pa;
double timer_rs_mv;
double timer_rs_un;
double timer_plot;
long long total_blocks;
num_sz nb_min;
num_sz nb_max;
int nrrs;
int nrs;
int nps;
int nlbs;
int num_refined;
int num_reformed;
int num_moved_all;
int num_moved_lb;
int num_moved_rs;
int num_moved_reduce;
int num_moved_coarsen;
int num_comm_x;
int num_comm_y;
int num_comm_z;
int num_comm_tot;
int num_comm_uniq;
int num_comm_x_min;
int num_comm_y_min;
int num_comm_z_min;
int num_comm_t_min;
int num_comm_u_min;
int num_comm_x_max;
int num_comm_y_max;
int num_comm_z_max;
int num_comm_t_max;
int num_comm_u_max;
int counter_halo_recv[3];
int counter_halo_send[3];
double size_mesg_recv[3];
double size_mesg_send[3];
int counter_face_recv[3];
int counter_face_send[3];
int counter_bc[3];
int counter_same[3];
int counter_diff[3];
int counter_malloc;
double size_malloc;
int counter_malloc_init;
double size_malloc_init;
int total_red;
