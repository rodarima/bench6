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

extern double average[146];
extern double stddev[143];
extern double minimum[143];
extern double maximum[143];

extern double timer_all;

extern double timer_comm_all;
extern double timer_comm_dir[3];
extern double timer_comm_recv[3];
extern double timer_comm_pack[3];
extern double timer_comm_send[3];
extern double timer_comm_same[3];
extern double timer_comm_diff[3];
extern double timer_comm_bc[3];
extern double timer_comm_wait[3];
extern double timer_comm_unpack[3];

extern double timer_calc_all;

extern double timer_cs_all;
extern double timer_cs_red;
extern double timer_cs_calc;

extern double timer_refine_all;
extern double timer_refine_co;
extern double timer_refine_mr;
extern double timer_refine_cc;
extern double timer_refine_sb;
extern double timer_refine_c1;
extern double timer_refine_c2;
extern double timer_refine_sy;
extern double timer_cb_all;
extern double timer_cb_cb;
extern double timer_cb_pa;
extern double timer_cb_mv;
extern double timer_cb_un;
extern double timer_lb_all;
extern double timer_lb_sort;
extern double timer_lb_pa;
extern double timer_lb_mv;
extern double timer_lb_un;
extern double timer_lb_misc;
extern double timer_lb_mb;
extern double timer_lb_ma;
extern double timer_rs_all;
extern double timer_rs_ca;
extern double timer_rs_pa;
extern double timer_rs_mv;
extern double timer_rs_un;

extern double timer_plot;

extern long long total_blocks;
extern num_sz nb_min;
extern num_sz nb_max;
extern int nrrs;
extern int nrs;
extern int nps;
extern int nlbs;
extern int num_refined;
extern int num_reformed;
extern int num_moved_all;
extern int num_moved_lb;
extern int num_moved_rs;
extern int num_moved_reduce;
extern int num_moved_coarsen;
extern int num_comm_x;
extern int num_comm_y;
extern int num_comm_z;
extern int num_comm_tot;
extern int num_comm_uniq;
extern int num_comm_x_min;
extern int num_comm_y_min;
extern int num_comm_z_min;
extern int num_comm_t_min;
extern int num_comm_u_min;
extern int num_comm_x_max;
extern int num_comm_y_max;
extern int num_comm_z_max;
extern int num_comm_t_max;
extern int num_comm_u_max;
extern int counter_halo_recv[3];
extern int counter_halo_send[3];
extern double size_mesg_recv[3];
extern double size_mesg_send[3];
extern int counter_face_recv[3];
extern int counter_face_send[3];
extern int counter_bc[3];
extern int counter_same[3];
extern int counter_diff[3];
extern int counter_malloc;
extern double size_malloc;
extern int counter_malloc_init;
extern double size_malloc_init;
extern int total_red;
