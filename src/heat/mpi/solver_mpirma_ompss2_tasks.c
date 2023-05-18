#include <mpi.h>

#include "utils.h"
#include "common/heat.h"


static int serial;

static inline void send(const double *data, int nelems, int dst, MPI_Win win)
{
	MPI_Win_fence(0, win);

	MPI_Put(data, nelems, MPI_DOUBLE, dst, 0, nelems, MPI_DOUBLE, win);

	MPI_Win_fence(0, win);
}

static inline void recv(MPI_Win win)
{
	MPI_Win_fence(0, win);

	MPI_Win_fence(0, win);
}

static inline void gaussSeidelSolver(
	int64_t rows, int64_t cols, int rbs, int cbs, int nrb, int ncb,
	int nwins, int ncbxw, MPI_Win upperWins[nwins], MPI_Win lowerWins[nwins],
	double M[rows][cols], char reps[nrb][ncb]
) {
	if (rank != 0) {
		for (int w = 0; w < nwins; ++w)
			#pragma oss task label("send upper inner") in({reps[1][1+C], C=w*ncbxw;ncbxw}) inout(serial)
			send(&M[1][1+w*ncbxw*cbs], ncbxw*cbs, 0, upperWins[w]);

		for (int w = 0; w < nwins; ++w)
			#pragma oss task label("recv upper outer") out({reps[0][1+C], C=w*ncbxw;ncbxw}) inout(serial)
			recv(upperWins[w]);
	}

	if (rank != nranks-1) {
		for (int w = 0; w < nwins; ++w)
			#pragma oss task label("recv lower outer") out({reps[nrb-1][1+C], C=w*ncbxw;ncbxw}) inout(serial)
			recv(lowerWins[w]);
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
		for (int w = 0; w < nwins; ++w)
			#pragma oss task label("send lower inner") in({reps[nrb-2][1+C], C=w*ncbxw;ncbxw}) inout(serial)
			send(&M[rows-2][1+w*ncbxw*cbs], ncbxw*cbs, 1, lowerWins[w]);
	}
}

double solve(HeatConfiguration *conf, int64_t rows, int64_t cols, int timesteps, void *extraData)
{
	double (*matrix)[cols] = (double (*)[cols]) conf->matrix;
	const int rbs = conf->rbs;
	const int cbs = conf->cbs;

	const int nrb = (rows-2)/rbs+2;
	const int ncb = (cols-2)/cbs+2;
	char representatives[nrb][ncb];

	const MPIRMAInfo *info = (MPIRMAInfo *) extraData;
	const int nwins = info->nwins;
	const int ncbxw = (ncb-2)/nwins;

	for (int t = 0; t < timesteps; ++t) {
		gaussSeidelSolver(rows, cols, rbs, cbs, nrb, ncb,
			nwins, ncbxw, info->upperWins, info->lowerWins,
			matrix, representatives);
	}
	#pragma oss taskwait

	MPI_Barrier(MPI_COMM_WORLD);

	return IGNORE_RESIDUAL;
}
