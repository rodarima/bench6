#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>

#include "utils.h"
#include "common/heat.h"

typedef struct {
	MPI_Request send;
	MPI_Request recv;
} HaloRequests;

static inline void put(const double *data, int cb, int cbs, int dst, MPI_Win win)
{
	MPI_Put(data, cbs, MPI_DOUBLE, dst, cb*cbs, cbs, MPI_DOUBLE, win);
}

static inline void irecv(int src, int tag, HaloRequests *reqs, MPI_Comm comm)
{
	MPI_Irecv(NULL, 0, MPI_BYTE, src, tag, comm, &reqs->recv);
}

static inline void inotify(int cb, int dst, HaloRequests *reqs, MPI_Comm comm, MPI_Win win)
{
	MPI_Win_flush(dst, win);
	MPI_Isend(NULL, 0, MPI_BYTE, dst, cb, comm, &reqs->send);
}

static inline void waitHaloSendRecv(int nhalos, HaloRequests *reqs)
{
	MPI_Waitall(2*nhalos, (MPI_Request *)reqs, MPI_STATUSES_IGNORE);
}

double solve(HeatConfiguration *conf, int64_t rows, int64_t cols, int timesteps, void *extraData)
{
	double (*M)[cols] = (double (*)[cols]) conf->matrix;
	const int rbs = conf->rbs;
	const int cbs = conf->cbs;
	const int nrb = (rows-2)/rbs+2;
	const int ncb = (cols-2)/cbs+2;

	const MPIRMAInfo *info = (MPIRMAInfo *) extraData;
	const int nwins = info->nwins;
	const int ncbxw = (ncb-2)/nwins;

	HaloRequests upperHaloReqs[nwins];
	HaloRequests lowerHaloReqs[nwins];
	for (int win = 0; win < nwins; ++win) {
		upperHaloReqs[win].send = MPI_REQUEST_NULL;
		upperHaloReqs[win].recv = MPI_REQUEST_NULL;
		lowerHaloReqs[win].send = MPI_REQUEST_NULL;
		lowerHaloReqs[win].recv = MPI_REQUEST_NULL;
	}

	MPI_Win_lock_all(0, info->upperWins[0]);
	MPI_Win_lock_all(0, info->lowerWins[0]);

	if (rank > 0)
		for (int win = 0; win < nwins; ++win)
			irecv(0, win, &upperHaloReqs[win], info->upperComm);

	for (int t = 0; t < timesteps; ++t) {
		for (int R = 1; R < nrb-1; ++R) {
			for (int C = 1; C < ncb-1; ++C) {
				// Compute the logical window of this block
				int win = (C-1)/ncbxw;
				int posinwin = (C-1)%ncbxw;

				// Wait for upper halo
				if (posinwin == 0 && R == 1 && rank > 0)
					waitHaloSendRecv(1, &upperHaloReqs[win]);

				// First row of the block
				const int fr = (R-1)*rbs+1;
				for (int c = (C-1)*cbs+1; c <= C*cbs; ++c) {
					M[fr][c] = 0.25*(M[fr-1][c] + M[fr+1][c] + M[fr][c-1] + M[fr][c+1]);
				}

				// Start sending upper compute row and receiving next upper halo
				if (fr == 1 && rank > 0) {
					put(&M[1][(C-1)*cbs+1], C-1, cbs, 0, info->upperWins[0]);
					if (posinwin == ncbxw-1)
						inotify(win, 0, &upperHaloReqs[win], info->upperComm, info->upperWins[0]);
					if (posinwin == 0 && t < timesteps-1)
						irecv(0, win, &upperHaloReqs[win], info->upperComm);
				}

				computeBlock(rows, cols, (R-1)*rbs+2, R*rbs-1, (C-1)*cbs+1, C*cbs, M);

				// Wait for lower halo
				const int lr = R*rbs;
				if (posinwin == 0 && lr == rows-2 && rank < nranks-1)
					waitHaloSendRecv(1, &lowerHaloReqs[win]);

				// Last row of the block
				for (int c = (C-1)*cbs+1; c <= C*cbs; ++c) {
					M[lr][c] = 0.25*(M[lr-1][c] + M[lr+1][c] + M[lr][c-1] + M[lr][c+1]);
				}

				// Start sending lower compute row and receiving next lower halo
				if (R == nrb-2 && rank < nranks-1) {
					put(&M[rows-2][(C-1)*cbs+1], C-1, cbs, 1, info->lowerWins[0]);
					if (posinwin == ncbxw-1)
						inotify(win, 1, &lowerHaloReqs[win], info->lowerComm, info->lowerWins[0]);
					if (posinwin == 0)
						irecv(1, win, &lowerHaloReqs[win], info->lowerComm);
				}
			}
		}
	}
	waitHaloSendRecv(nwins, upperHaloReqs);
	waitHaloSendRecv(nwins, lowerHaloReqs);

	MPI_Win_unlock_all(info->upperWins[0]);
	MPI_Win_unlock_all(info->lowerWins[0]);

	MPI_Barrier(MPI_COMM_WORLD);

	return IGNORE_RESIDUAL;
}
