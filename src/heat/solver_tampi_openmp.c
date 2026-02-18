#include <mpi.h>
#include <TAMPI.h>

#include "utils_mpi.h"
#include "heat.h"

int
mpi_level(void)
{
	return MPI_THREAD_MULTIPLE;
}

const char *
summary(void)
{
	return "Parallel version using OpenMP tasks and TAMPI";
}

static inline void send(const double *data, int nelems, int dst, int tag)
{
#if !defined(TAMPI_VERSION_MAJOR) || (TAMPI_VERSION_MAJOR == 3)
	MPI_Request request;
	MPI_Isend(data, nelems, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &request);
	TAMPI_Iwait(&request, MPI_STATUS_IGNORE);
#elif (TAMPI_VERSION_MAJOR == 4)
	TAMPI_Isend(data, nelems, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
#else
	#error "TAMPI version not supported for this benchmark"
#endif
}

static inline void recv(double *data, int nelems, int src, int tag)
{
#if !defined(TAMPI_VERSION_MAJOR) || (TAMPI_VERSION_MAJOR == 3)
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
		for (int C = 1; C < ncb-1; ++C) {
			#pragma omp task depend(in: reps[1][C]) priority(1) label("send upper inner")
			send(&M[1][(C-1)*cbs+1], cbs, rank-1, C);
		}
		for (int C = 1; C < ncb-1; ++C) {
			#pragma omp task depend(out: reps[0][C]) priority(1) label("recv upper outer")
			recv(&M[0][(C-1)*cbs+1], cbs, rank-1, C);
		}
    }

	if (rank != nranks-1) {
		for (int C = 1; C < ncb-1; ++C) {
			#pragma omp task depend(out: reps[nrb-1][C]) priority(1) label("recv lower outer")
			recv(&M[rows-1][(C-1)*cbs+1], cbs, rank+1, C);
		}
	}

	for (int R = 1; R < nrb-1; ++R) {
		for (int C = 1; C < ncb-1; ++C) {
			#pragma omp task \
				depend(in: reps[R+1][C]) \
				depend(in: reps[R-1][C]) \
				depend(in: reps[R][C-1]) \
				depend(in: reps[R][C+1]) \
				depend(inout: reps[R][C]) \
				label("block computation")
			computeBlock(rows, cols, (R-1)*rbs+1, R*rbs, (C-1)*cbs+1, C*cbs, M);
		}
	}

	if (rank != nranks-1) {
		for (int C = 1; C < ncb-1; ++C) {
			#pragma omp task depend(in: reps[nrb-2][C]) priority(1) label("send lower inner")
			send(&M[rows-2][(C-1)*cbs+1], cbs, rank+1, C);
		}
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

	#pragma omp parallel default(none) shared(representatives, matrix, rows, cols, rbs, cbs, nrb, ncb) firstprivate(timesteps)
	{
		#pragma omp single nowait
		{
			for (int t = 0; t < timesteps; ++t) {
				gaussSeidelSolver(rows, cols, rbs, cbs, nrb, ncb, matrix, representatives);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	return IGNORE_RESIDUAL;
}
