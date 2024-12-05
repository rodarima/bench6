#include <mpi.h>
#include <TAMPI.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "common/matmul.h"

int rank, nranks;


int main(int argc, char **argv)
{
	// Initialize MPI and TAMPI
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_TASK_MULTIPLE, &provided);
	if (provided != MPI_TASK_MULTIPLE)
		matmul_fail("MPI_TASK_MULTIPLE not supported");

	MPI_Info info;
	MPI_Info_create(&info);
	MPI_Info_set(info, "mpi_assert_no_any_tag", "true");
	MPI_Info_set(info, "mpi_assert_no_any_source", "true");
	MPI_Info_set(info, "mpi_assert_exact_length", "true");
	MPI_Comm_set_info(MPI_COMM_WORLD, info);
	MPI_Info_free(&info);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);

	matmul_conf_t conf;
	if (!rank)
		matmul_getconf(argc, argv, &conf);

	MPI_Bcast(&conf, sizeof(conf), MPI_BYTE, 0, MPI_COMM_WORLD);

	const size_t n = conf.n;
	const size_t m_per_rank = conf.m / nranks;
	const size_t ts = conf.ts;

	if (conf.n % nranks)
		matmul_fail("The matrix size must be divisible by the number of processes");
	if (conf.n % conf.ts != 0 || m_per_rank % conf.ts != 0)
		matmul_fail("The matrix size must be divisible by the tile size");

	matmul_t matmul;
	matmul_setup(&conf, &matmul, m_per_rank);

	MPI_Barrier(MPI_COMM_WORLD);

	if (conf.warmup) {
		matmul_solve(n, m_per_rank, ts, &matmul, conf.warmup);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	// Execute matmul solver
	double start = matmul_gettime();

	matmul_solve(n, m_per_rank, ts, &matmul, conf.timesteps);
	#pragma oss taskwait

	MPI_Barrier(MPI_COMM_WORLD);

	double end = matmul_gettime();

	double expected = (double)(conf.timesteps + conf.warmup) * matmul.alpha * matmul.beta * n;

	// Check results
	int errors = matmul_check(n, m_per_rank, ts, matmul.C, expected);

	MPI_Allreduce(MPI_IN_PLACE, &errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if (!rank && errors > 0)
		fprintf(stderr, "Results checked: WRONG, mismatches: %d\n", errors);
	else if (!rank) {
		fprintf(stdout, "Results checked: CORRECT\n");
		fprintf(stdout, "n: %zu, ts: %zu, timesteps: %zu, time(s): %f\n", n, ts, conf.timesteps, end-start);
	}

	// Free memory
	free(matmul.A);
	free(matmul.remote1);
	free(matmul.remote2);
	free(matmul.B);
	free(matmul.C);

	// Finalize MPI and TAMPI
	MPI_Finalize();

	return 0;
}
