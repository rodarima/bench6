#include <mpi.h>

#include "utils_mpi.h"
#include "heat.h"

int
mpi_level(void)
{
	return MPI_THREAD_SERIALIZED;
}

const char *
summary(void)
{
	return "Parallel version using MPI + OmpSs-2 following a fork-join\n"
		"parallelization.";
}

static inline void send(const double *data, int nelems, int dst, int tag)
{
	MPI_Send(data, nelems, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
}

static inline void recv(double *data, int nelems, int src, int tag)
{
	MPI_Recv(data, nelems, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

static inline void gaussSeidelSolver(int64_t rows, int64_t cols, int rbs, int cbs, double M[rows][cols])
{
	const int nrb = (rows-2)/rbs+2;
	const int ncb = (cols-2)/cbs+2;
	char reps[nrb][ncb];

	if (rank != 0) {
		for (int C = 1; C < cols-1; C += cbs)
			send(&M[1][C], cbs, rank-1, C);
		for (int C = 1; C < cols-1; C += cbs)
			recv(&M[0][C], cbs, rank-1, C);
    }

	if (rank != nranks-1) {
		for (int C = 1; C < cols-1; C += cbs)
			recv(&M[rows-1][C], cbs, rank+1, C);
	}

	for (int R = 1; R < nrb-1; ++R) {
		for (int C = 1; C < ncb-1; ++C) {
			#pragma oss task label("block computation") \
					in(reps[R-1][C]) in(reps[R+1][C]) \
					in(reps[R][C-1]) in(reps[R][C+1]) \
					inout(reps[R][C])
			computeBlock(rows, cols, (R-1)*rbs+1, R*rbs, (C-1)*cbs+1, C*cbs, M);
		}
	}
	#pragma oss taskwait

	if (rank != nranks-1) {
		for (int C = 1; C < cols-1; C += cbs)
			send(&M[rows-2][C], cbs, rank+1, C);
	}
}

double solve(HeatConfiguration *conf, int64_t rows, int64_t cols, int timesteps, void *extraData)
{
	(void) extraData;
	double (*matrix)[cols] = (double (*)[cols]) conf->matrix;

	for (int t = 0; t < timesteps; ++t) {
		gaussSeidelSolver(rows, cols, conf->rbs, conf->cbs, matrix);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	return IGNORE_RESIDUAL;
}
