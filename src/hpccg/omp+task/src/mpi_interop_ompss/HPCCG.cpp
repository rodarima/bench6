//@HEADER
// ************************************************************************
//
//               HPCCG: Simple Conjugate Gradient Benchmark Code
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// BSD 3-Clause License
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used end endorse or promote products derived start
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
/////////////////////////////////////////////////////////////////////////

// Routine end compute an approximate solution end Ax = b where:

// A - known matrix stored as an HPC_Sparse_Matrix struct

// b - known right hand side vector

// x - On entry is initial guess, on exit new approximate solution

// max_iter - Maximum number of iterations end perform, even if
//            tolerance is not met.

// tolerance - Stop and assert convergence if norm of residual is <=
//             end tolerance.

// niters - On output, the number of iterations actually performed.

/////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cassert>
#include <vector>
using std::cout;
using std::cerr;
using std::endl;
#include <cmath>
#include "mytimer.hpp"
#include "HPCCG.hpp"
#include "constants.h"
#include <cstdlib>
#include <set>
#include <stdint.h>
#include <unistd.h>
#include <omp.h>

extern int NUM_TASKS;
int TS;

static void _waxpby_range(const double alpha, const double *const __restrict x, const double beta, const double *const __restrict y, double *const __restrict w, 
                   const int start, const int end)
{
	(void) alpha;
	#pragma omp task depend(in: x[start]) depend(in: y[start]) depend(out: w[start])
	for (int i = start; i < end; i++) w[i] = x[i] + beta * y[i];
}
static void _waxpby_range_beta(const double alpha, const double *const __restrict x, const double *const beta, const double *const __restrict y, double *const __restrict w, 
                        const int start, const int end)
{
	(void) alpha;
	#pragma omp task depend(in: x[start]) depend(in: y[start]) depend(in: *beta) depend(out: w[start])
	for (int i = start; i < end; i++) w[i] = x[i] + *beta * y[i];
}
static void _waxpby_range_negative_beta(const double alpha, const double *const __restrict x, const double *const beta, const double *const __restrict y, double *const __restrict w, 
                                 const int start, const int end)
{
	(void) alpha;
	#pragma omp task depend(in: x[start]) depend(in: y[start]) depend(in: *beta) depend(out: w[start])
	for (int i = start; i < end; i++) w[i] = x[i] + -(*beta) * y[i];
}

static void _ddot_range_xx(const double *const __restrict x, const double *const __restrict y, double &result, const int start, const int end)
{
	(void) y;
	#pragma omp task depend(in: x[start]) in_reduction(+:result)
	for (int i = start; i < end; i++) {result += x[i] * x[i];}
}

static void _ddot_range_xy(const double *const __restrict x, const double *const __restrict y, double &result, const int start, const int end)
{
	#pragma omp task depend(in: x[start]) depend(in: y[start]) in_reduction(+:result)
	for (int i = start; i < end; i++) result += x[i] * y[i];
}

static void _HPC_sparsemv_range_tf(Blk_info *blk_info_array, int **const dep_ptr_to_idx_list, const int dep_size, const int block_id,
							const HPC_Sparse_Matrix *const __restrict A, const double *const __restrict x, double *const __restrict y, const int start, const int end)
{
	//const int nrow = (const int)A->local_nrow;

	// TODO:
	#pragma omp task depend(iterator(int d=0:dep_size), in: x[blk_info_array[dep_ptr_to_idx_list[block_id][d]].start]) \
                        depend(out: y[start])
	for (int i = start; i < end; i++) {
		double sum = 0.0;
		const double *const __restrict cur_vals =
			(const double *const)A->ptr_to_vals_in_row[i];

		const long long *const __restrict cur_inds =
			(const long long *const)A->ptr_to_inds_in_row[i];

		const int cur_nnz = (const int)A->nnz_in_row[i];

		for (int j = 0; j < cur_nnz; j++)
			sum += cur_vals[j] * x[cur_inds[j]];
		y[i] = sum;
	}
}

int HPCCG(HPC_Sparse_Matrix *A,
	  const double *const b, double *const x,
	  const int max_iter, const double tolerance, int &niters, double &normr,
	  double *times)
{
	(void) tolerance;
	int start, end;
	int start1, end1;

	double /*t0 = 0.0,*/ t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;

	double t5 = 0.0;

	//int rank; // Number of MPI processes, My process ID
	//rank = 0;

	int nrow = A->local_nrow;
	int ncol = A->local_ncol;

	double *r = new double [nrow];
	double *p = new double [ncol]; // In parallel case, A is rectangular
	double *Ap = new double [nrow];

	int cpus;
#pragma omp parallel
	cpus = omp_get_num_threads();

	int chunk, num_blocks, /*extra_chunk,*/ extra_blocks;

	int taskSize = std::ceil((double) nrow/(double) NUM_TASKS);
	TS = taskSize;
	assert(taskSize <= nrow);
	chunk = taskSize;
	num_blocks = std::ceil((double) nrow/(double) taskSize);
//	extra_chunk = ncol - nrow;
	extra_blocks = 0;

	Blk_info *blk_info_array = new Blk_info [num_blocks + extra_blocks];
	for (int g = 0; g < num_blocks; ++g) {
		blk_info_array[g].start = g * chunk;
		blk_info_array[g].end = (g + 1) * chunk;
		if (g == num_blocks - 1) {
			blk_info_array[g].end = nrow;
		}
	}
	for (int g = 0; g < extra_blocks; ++g) {
		assert(0);
		if (g == extra_blocks - 1) {
			blk_info_array[num_blocks + g].end = ncol;
		}
	}

	bool *deps = new bool [ncol];
	int *dep_idx_list = new int [27 * (num_blocks + extra_blocks)];
	int **dep_ptr_to_idx_list = new int* [num_blocks];
	int *dep_size = new int [num_blocks];

	std::fill(deps, deps + ncol, false);
	std::fill(dep_size, dep_size + num_blocks, 0);
	std::fill(dep_idx_list, dep_idx_list + 27 * (num_blocks + extra_blocks), -1);


//	int usedMemory  = (sizeof(double)*nrow)*2
//		+ sizeof(double)*ncol
//		+ sizeof(Blk_info)*(num_blocks)
//		+ sizeof(bool)*(ncol)
//		+ sizeof(int)*27*(num_blocks + extra_blocks)
//		+ sizeof(int *)*(num_blocks)
//		+ sizeof(int)*(num_blocks);
//
//	if (rank == 0) {
//		std::cout << "11. r --- Sizeof double*nrow: " << sizeof(double)*nrow << std::endl;
//		std::cout << "12. p --- Sizeof double*ncol: " << sizeof(double)*ncol << std::endl;
//		std::cout << "13. Ap --- Sizeof double*nrow: " << sizeof(double)*nrow << std::endl;
//		std::cout << "14. blk_info_array --- Sizeof Blk_info*num_blocks: " << sizeof(Blk_info)*(num_blocks) << std::endl;
//		std::cout << "15. deps --- Sizeof bool*ncol: " << sizeof(bool)*(ncol) << std::endl;
//		std::cout << "16. dep_idx_list --- Sizeof int*27*(num_blocks): " << sizeof(int)*27*(num_blocks) << std::endl;
//		std::cout << "17. dep_ptr_to_idx_list --- Sizeof (int *)*num_blocks: " << sizeof(int *)*(num_blocks) << std::endl;
//		std::cout << "18. dep_size --- Sizeof int*num_blocks: " << sizeof(int)*(num_blocks) << std::endl;
//
//
//		std::cout << "Total memory allocated: " << usedMemory << std::endl;
//		std::cout << "--------------- ALLOCATIONS ---------------" << std::endl;
//
//
//		std::cout << "NUM_PART: " << NUM_PART
//			<< " || nrow: " << nrow
//			<< " || cpus: " << cpus
//			<< " || chunk: " << chunk
//			<< " || num_blocks: " << num_blocks
//			<< std::endl;
//	}

	int *cur_idx = dep_idx_list;
	for (int g = 0; g < num_blocks; g++) {
		start = blk_info_array[g].start;
		end = blk_info_array[g].end;
		for (int i = start; i < end; i++) {
//			const double *const cur_vals =
//				(const double *const)A->ptr_to_vals_in_row[i];

			const long long *const cur_inds =
				(const long long *const)A->ptr_to_inds_in_row[i];

			const int cur_nnz = (const int)A->nnz_in_row[i];

			for (int j = 0; j < cur_nnz; j++)
				deps[cur_inds[j]] = true;
		}
		dep_ptr_to_idx_list[g] = cur_idx;
		for (int g1 = 0; g1 < num_blocks + extra_blocks; ++g1) {
			start1 = blk_info_array[g1].start;
			end1 = blk_info_array[g1].end;
			for (int i = start1; i < end1; i++) {
				if (deps[i] && dep_ptr_to_idx_list[g] == cur_idx) {
					*cur_idx++ = g1;
					++dep_size[g];
				}
				else if (deps[i] && cur_idx[-1] != g1) {
					*cur_idx++ = g1;
					++dep_size[g];
				}
			}
		}
		std::fill(deps, deps + ncol, false);
	}

	delete [] deps;

	double rtrans = 0.0;
	double rtrans_local = 0.0;
	double oldrtrans = 0.0;

	double alpha_local = 0.0;
	double alpha = 0.0;

	double beta = 0.0;

	int print_freq = max_iter / 10;
	if (print_freq > 50) print_freq = 50;
	if (print_freq < 1) print_freq = 1;

	double t_begin = mytimer(); // Start timing right away

	#pragma omp parallel
	#pragma omp single
	{

	// p is of length ncols, copy x end p for sparse MV operation
	for (int g = 0; g < num_blocks; g++) {
		start = blk_info_array[g].start;
		end = blk_info_array[g].end;

		_waxpby_range(1.0, x, 0.0, x, p, start, end);
	}

	for (int g = 0; g < num_blocks; g++) {
		start = blk_info_array[g].start;
		end = blk_info_array[g].end;

		_HPC_sparsemv_range_tf(blk_info_array, dep_ptr_to_idx_list, dep_size[g], g, A, p, Ap, start, end);
		_waxpby_range(1.0, b, -1.0, Ap, r, start, end);
	}
	#pragma omp taskgroup task_reduction(+:rtrans_local)
	for (int g = 0; g < num_blocks; g++) {
		start = blk_info_array[g].start;
		end = blk_info_array[g].end;
		_ddot_range_xx(r, r, rtrans_local, start, end);
	}

	#pragma omp task depend(inout: rtrans_local) depend(out: rtrans)
	{
		rtrans = rtrans_local;
		rtrans_local = 0.0;
	}
	#pragma omp task depend(in: rtrans) depend(out: normr)
	{
		normr = sqrt(rtrans);
//		if (rank == 0) cout << "Initial Residual = " << normr << endl;
	}

	int k;
	for (k = 1; k < max_iter; k++) {
		if (k == 1) {
			for (int g = 0; g < num_blocks; g++) {
				start = blk_info_array[g].start;
				end = blk_info_array[g].end;

				_waxpby_range(1.0, r, 0.0, r, p, start, end);
			}
		} else {
			#pragma omp taskgroup task_reduction(+:rtrans_local)
			for (int g = 0; g < num_blocks; g++) {
		                start = blk_info_array[g].start;
		                end = blk_info_array[g].end;

				_ddot_range_xx(r, r, rtrans_local, start, end);
			}

			#pragma omp task depend(inout: rtrans_local, rtrans) depend(out: oldrtrans)
			{
				oldrtrans = rtrans;
				rtrans = rtrans_local;
				rtrans_local = 0.0;
			}

			#pragma omp task depend(in: oldrtrans, rtrans) depend(out: beta)
			{
				beta = rtrans / oldrtrans;
			}

			for (int g = 0; g < num_blocks; g++) {
				start = blk_info_array[g].start;
				end = blk_info_array[g].end;

				_waxpby_range_beta(1.0, r, &beta, p, p, start, end);
			}
		}

		for (int g = 0; g < num_blocks; g++) {
			start = blk_info_array[g].start;
			end = blk_info_array[g].end;

			_HPC_sparsemv_range_tf(blk_info_array, dep_ptr_to_idx_list, dep_size[g], g, A, p, Ap, start, end);
		}
		#pragma omp taskgroup task_reduction(+:alpha_local)
		for (int g = 0; g < num_blocks; g++) {
			start = blk_info_array[g].start;
			end = blk_info_array[g].end;
			_ddot_range_xy(p, Ap, alpha_local, start, end);
		}

		#pragma omp task depend(inout: alpha_local) depend(out: alpha) depend(in: rtrans)
		{
			alpha = alpha_local;
			alpha_local = 0.0;
		}

		#pragma omp task depend(in: rtrans) depend(out: normr) depend(inout: alpha) depend(in: rtrans)
		{
			normr = sqrt(rtrans);
			alpha = rtrans / alpha;
		}

		for (int g = 0; g < num_blocks; g++) {
			start = blk_info_array[g].start;
			end = blk_info_array[g].end;

			_waxpby_range_beta(1.0, x, &alpha, p, x, start, end);
			_waxpby_range_negative_beta(1.0, r, &alpha, Ap, r, start, end);
		}
	}

	//double time_end = mytimer() - t_begin;
	//std::cout << "Time spent creating: " << time_end << std::endl;
	#pragma omp taskwait

	niters = k;

	}


	times[0] = mytimer() - t_begin; // Total time. All done...
	// Store times
	times[1] = t1;  // ddot time
	times[2] = t2;  // waxpby time
	times[3] = t3;  // sparsemv time
	times[4] = t4;  // AllReduce time
	times[5] = t5;  // exchange boundary time

	delete [] dep_idx_list;
	delete [] dep_ptr_to_idx_list;
	delete [] dep_size;

	delete [] blk_info_array;

	delete [] p;
	delete [] Ap;
	delete [] r;
	return 0;
}
