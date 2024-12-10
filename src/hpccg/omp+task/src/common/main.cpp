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
//   contributors may be used to endorse or promote products derived from
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

// Main routine of a program that reads a sparse matrix, right side
// vector, solution vector and initial guess from a file  in HPC
// format.  This program then calls the HPCCG conjugate gradient
// solver to solve the problem, and then prints results.

// Calling sequence:

// test_HPCCG linear_system_file

// Routines called:

// read_HPC_row - Reads in linear system

// mytimer - Timing routine (compile with -DWALL to get wall clock
//           times

// HPCCG - CG Solver

// compute_residual - Compares HPCCG solution to known solution.

#include <iostream>
#include <sstream>
using std::cout;
using std::cerr;
using std::endl;
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cassert>
#include <string>
#include <cmath>
#include <omp.h>
#include "generate_matrix.hpp"
#include "read_HPC_row.hpp"
#include "mytimer.hpp"
#include "HPC_sparsemv.hpp"
#include "compute_residual.hpp"
#include "HPCCG.hpp"
#include "HPC_Sparse_Matrix.hpp"
#include "dump_matlab_matrix.hpp"

#include "YAML_Element.hpp"
#include "YAML_Doc.hpp"

#undef DEBUG

int NUM_TASKS;

int main(int argc, char *argv[])
{
	HPC_Sparse_Matrix *A;
	double *x, *b, *xexact;
	int ierr = 0;
	double times[7];
	int nx = 0, ny = 0, nz = 0;
	int max_iter = 150;
	int num_tasks;

	int size = 1; // Serial case (not using MPI)
	int rank = 0;

#ifdef DEBUG
	if (rank == 0) {
		int junk = 0;
		cout << "Press enter to continue" << endl;
		cin >> junk;
	}

	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if (argc != 2 && argc != 6) {
		if (rank == 0) {
			cerr << "Usage:" << endl
				<< "Mode 1: " << argv[0] << " nx ny nz max_it num_tasks" << endl
				<< "     where nx, ny and nz are the local sub-block dimensions, max_it the" << endl
				<< "     maximum number of iterations of the algorithm, " << endl
				<< "     num_tasks the number of tasks per kernel (only for tasking/tf versions)" << endl
				<< "Mode 2: " << argv[0] << " HPC_data_file " << endl
				<< "     where HPC_data_file is a globally accessible file containing matrix data." << endl;
		}
		exit(1);
	}

	if (argc == 6) {
		nx = atoi(argv[1]);
		ny = atoi(argv[2]);
		nz = atoi(argv[3]);
		max_iter = atoi(argv[4]);
		num_tasks = atoi(argv[5]);
		// TODO: remove NUM_TASKS
		NUM_TASKS = num_tasks;
		generate_matrix(nx, ny, nz, &A, &x, &b, &xexact);
	} else {
		read_HPC_row(argv[1], &A, &x, &b, &xexact);
	}


	bool dump_matrix = false;
	if (dump_matrix && size <= 4) dump_matlab_matrix(A, rank);

	//double t1 = mytimer(); // Initialize it (if needed)
	int niters = 0;
	double normr = 0.0;
	double tolerance = 0.0; // Set tolerance to zero to make all runs do max_iter iterations

	// ierr = HPCCG(A, b, x, max_iter, tolerance, niters, normr, times, A->local_nrow/num_tasks);
	ierr = HPCCG(A, b, x, max_iter, tolerance, niters, normr, times);

	if (ierr) cerr << "Error in call to CG: " << ierr << ".\n" << endl;

	// initialize YAML doc
	if (rank == 0) { // Only PE 0 needs to compute and report timing results
		double fniters = niters;
		double fnrow = A->total_nrow;
		double fnnz = A->total_nnz;
		double fnops_ddot = fniters * 4 * fnrow;
		double fnops_waxpby = fniters * 6 * fnrow;
		double fnops_sparsemv = fniters * 2 * fnnz;
		double fnops = fnops_ddot + fnops_waxpby + fnops_sparsemv;

		YAML_Doc doc("hpccg", "1.0");
		doc.add("Parallelism", "");

		doc.get("Parallelism")->add("MPI not enabled", "");

		int nthreads = 1;
#pragma omp parallel
		nthreads = omp_get_num_threads();
		doc.get("Parallelism")->add("Number of OpenMP threads", nthreads);

		doc.add("Dimensions", "");
		doc.get("Dimensions")->add("nx", (int) nx);
		doc.get("Dimensions")->add("ny", (int) ny);
		doc.get("Dimensions")->add("nz", (int) nz);

		doc.add("Number of iterations", (int) niters);
		doc.add("Final residual", normr);
		doc.add("#********** Performance Summary (times in sec) ***********", "");

		doc.add("Time Summary", "");
		doc.get("Time Summary")->add("Total   ", times[0]);

		doc.add("FLOPS Summary", "");
		doc.get("FLOPS Summary")->add("Total   ", fnops);

		doc.add("MFLOPS Summary", "");
		doc.get("MFLOPS Summary")->add("Total   ", fnops / times[0] / 1.0E6);

//		if (rank == 0) { // only PE 0 needs to compute and report timing results
//			std::string yaml = doc.generateYAML();
//			cout << yaml;
//		}
		printf("%14e %d %d %d %d %d\n", times[0],
				nx, ny, nz, max_iter, num_tasks);
	}

	// Compute difference between known exact solution and computed solution
	// All processors are needed here.

	//double residual = 0;
	//  if ((ierr = compute_residual(A->local_nrow, x, xexact, &residual)))
	//  cerr << "Error in call to compute_residual: " << ierr << ".\n" << endl;

	// if (rank==0)
	//   cout << "Difference between computed and exact  = "
	//        << residual << ".\n" << endl;


	// Finish up
	delete_matrix(&A, &x, &b, &xexact, nx, ny, nz);
	return 0;
}
