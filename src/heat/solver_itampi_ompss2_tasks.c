#include <mpi.h>
#include <TAMPI.h>

#include "utils_mpi.h"
#include "heat.h"

const char *
summary(void)
{
	return "Parallel version using MPI + OmpSs-2 tasks + Non-blocking TAMPI";
}

int
mpi_level(void)
{
	return MPI_TASK_MULTIPLE;
}

static inline void send(const double *data, int nelems, int dst, int tag)
{
	#if (TAMPI_VERSION_MAJOR == 4)
	TAMPI_Isend(data, nelems, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
	#elif (TAMPI_VERSION_MAJOR == 3)
	MPI_Request request;
	MPI_Isend(data, nelems, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &request);
	TAMPI_Iwait(&request, MPI_STATUS_IGNORE);
	#else
	#error "TAMPI version not supported for this benchmark"
	#endif
}

static inline void recv(double *data, int nelems, int src, int tag)
{
#if (TAMPI_VERSION_MAJOR == 3)
	MPI_Request request;
	MPI_Irecv(data, nelems, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &request);
	TAMPI_Iwait(&request, MPI_STATUS_IGNORE);
#elif (TAMPI_VERSION_MAJOR == 4)
	TAMPI_Irecv(data, nelems, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#else
	#error "TAMPI version not supported for this benchmark"
#endif
}

static inline void gaussSeidelSolver(int64_t rows, int64_t cols, int rbs, int cbs, int nrb, int ncb, double M[rows][cols], char reps[nrb][ncb])
{
	if (rank != 0) {
		for (int C = 1; C < ncb-1; ++C)
			#pragma oss task label("send upper inner") in(reps[1][C])
			send(&M[1][(C-1)*cbs+1], cbs, rank-1, C);
		for (int C = 1; C < ncb-1; ++C)
			#pragma oss task label("recv upper outer") out(reps[0][C])
			recv(&M[0][(C-1)*cbs+1], cbs, rank-1, C);
    }

	if (rank != nranks-1) {
		for (int C = 1; C < ncb-1; ++C)
			#pragma oss task label("recv lower outer") out(reps[nrb-1][C])
			recv(&M[rows-1][(C-1)*cbs+1], cbs, rank+1, C);
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

	if (rank != nranks-1) {
		for (int C = 1; C < ncb-1; ++C)
			#pragma oss task label("send lower inner") in(reps[nrb-2][C])
			send(&M[rows-2][(C-1)*cbs+1], cbs, rank+1, C);
	}
}

double solve(HeatConfiguration *conf, int64_t rows, int64_t cols, int timesteps, void *extraData)
{
	(void) extraData;
	double (*matrix)[cols] = (double (*)[cols]) conf->matrix;
	const int rbs = conf->rbs;
	const int cbs = conf->cbs;

	const int nrb = (rows-2)/rbs+2;
	const int ncb = (cols-2)/cbs+2;
	char representatives[nrb][ncb];

	for (int t = 0; t < timesteps; ++t) {
		gaussSeidelSolver(rows, cols, rbs, cbs, nrb, ncb, matrix, representatives);
	}
	#pragma oss taskwait

	MPI_Barrier(MPI_COMM_WORLD);

	return IGNORE_RESIDUAL;
}
