#include <mpi.h>
#include <TAMPI.h>

#include "utils.h"
#include "heat.h"


static inline void fence(MPI_Win win)
{
	MPI_Request request;
	MPI_Win_ifence(0, win, &request);
	MPI_Wait(&request, MPI_STATUS_IGNORE);
}

static inline void ifence(MPI_Win win)
{
	MPI_Request request;
	MPI_Win_ifence(0, win, &request);
	TAMPI_Iwait(&request, MPI_STATUS_IGNORE);
}

static inline void send(const double *data, int nelems, int roffset, int dst, MPI_Win win)
{
	MPI_Put(data, nelems, MPI_DOUBLE, dst, roffset, nelems, MPI_DOUBLE, win);
}

static inline void recv(MPI_Win win)
{
	fence(win);

	ifence(win);
}

static inline void gaussSeidelSolver(
	int64_t rows, int64_t cols, int rbs, int cbs, int nrb, int ncb,
	int nwins, int ncbxw, MPI_Win upperWins[nwins], MPI_Win lowerWins[nwins],
	double M[rows][cols], char reps[nrb][ncb]
) {
	if (rank != 0) {
		for (int w = 0; w < nwins; ++w) {
			#pragma oss task label("ifence upper inner") in({reps[1][1+C], C=w*ncbxw;ncbxw}) inout(upperWins[w])
			ifence(upperWins[w]);

			for (int C = w*ncbxw+1; C < (w+1)*ncbxw+1; ++C)
				#pragma oss task label("send upper inner") in(reps[1][C]) concurrent(upperWins[w])
				send(&M[1][(C-1)*cbs+1], cbs, (C-w*ncbxw-1)*cbs, 0, upperWins[w]);

			#pragma oss task label("ifence upper inner") in({reps[1][1+C], C=w*ncbxw;ncbxw}) inout(upperWins[w])
			ifence(upperWins[w]);
		}

		for (int w = 0; w < nwins; ++w)
			#pragma oss task label("recv upper outer") out({reps[0][1+C], C=w*ncbxw;ncbxw}) inout(upperWins[w])
			recv(upperWins[w]);
	}

	if (rank != nranks-1) {
		for (int w = 0; w < nwins; ++w)
			#pragma oss task label("recv lower outer") out({reps[nrb-1][1+C], C=w*ncbxw;ncbxw}) inout(lowerWins[w])
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
		for (int w = 0; w < nwins; ++w) {
			#pragma oss task label("ifence lower inner") in({reps[nrb-2][1+C], C=w*ncbxw;ncbxw}) inout(lowerWins[w])
			ifence(lowerWins[w]);

			for (int C = w*ncbxw+1; C < (w+1)*ncbxw+1; ++C)
				#pragma oss task label("send lower inner") in(reps[nrb-2][C]) concurrent(lowerWins[w])
				send(&M[rows-2][(C-1)*cbs+1], cbs, (C-w*ncbxw-1)*cbs, 1, lowerWins[w]);

			#pragma oss task label(ifence lower inner) in({reps[nrb-2][1+C], C=w*ncbxw;ncbxw}) inout(lowerWins[w])
			ifence(lowerWins[w]);
		}
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
