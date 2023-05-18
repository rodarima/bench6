#include <mpi.h>

#include "utils.h"
#include "common/heat.h"

typedef struct {
	MPI_Request send;
	MPI_Request recv;
} HaloRequests;

static inline void isend(const double *data, int nelems, int dst, int tag, HaloRequests *reqs)
{
	MPI_Isend(data, nelems, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &reqs->send);
}

static inline void irecv(double *data, int nelems, int src, int tag, HaloRequests *reqs)
{
	MPI_Irecv(data, nelems, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &reqs->recv);
}

static inline void waitHaloSendRecv(int nhalos, HaloRequests *reqs)
{
	MPI_Waitall(2*nhalos, (MPI_Request *)reqs, MPI_STATUSES_IGNORE);
}

double solve(HeatConfiguration *conf, int64_t rows, int64_t cols, int timesteps, void *extraData)
{
	(void) extraData;
	double (*M)[cols] = (double (*)[cols]) conf->matrix;
	const int rbs = conf->rbs;
	const int cbs = conf->cbs;
	const int nrb = (rows-2)/rbs+2;
	const int ncb = (cols-2)/cbs+2;

	HaloRequests upperHaloReqs[ncb];
	HaloRequests lowerHaloReqs[ncb];
	for (int C = 0; C < ncb; ++C) {
		upperHaloReqs[C].send = MPI_REQUEST_NULL;
		upperHaloReqs[C].recv = MPI_REQUEST_NULL;
		lowerHaloReqs[C].send = MPI_REQUEST_NULL;
		lowerHaloReqs[C].recv = MPI_REQUEST_NULL;
	}

	if (rank > 0)
		for (int C = 1; C < ncb-1; ++C)
			irecv(&M[0][(C-1)*cbs+1], cbs, rank-1, C, &upperHaloReqs[C]);

	for (int t = 0; t < timesteps; ++t) {
		for (int R = 1; R < nrb-1; ++R) {
			for (int C = 1; C < ncb-1; ++C) {
				if (R == 1 && rank > 0)
					waitHaloSendRecv(1, &upperHaloReqs[C]);

				// First row of the block
				const int fr = (R-1)*rbs+1;
				for (int c = (C-1)*cbs+1; c <= C*cbs; ++c) {
					M[fr][c] = 0.25*(M[fr-1][c] + M[fr+1][c] + M[fr][c-1] + M[fr][c+1]);
				}

				if (fr == 1 && rank > 0) {
					isend(&M[1][(C-1)*cbs+1], cbs, rank-1, C, &upperHaloReqs[C]);
					if (t < timesteps-1)
						irecv(&M[0][(C-1)*cbs+1], cbs, rank-1, C, &upperHaloReqs[C]);
				}

				computeBlock(rows, cols, (R-1)*rbs+2, R*rbs-1, (C-1)*cbs+1, C*cbs, M);

				const int lr = R*rbs;
				if (lr == rows-2 && rank < nranks-1)
					waitHaloSendRecv(1, &lowerHaloReqs[C]);

				// Last row of the block
				for (int c = (C-1)*cbs+1; c <= C*cbs; ++c) {
					M[lr][c] = 0.25*(M[lr-1][c] + M[lr+1][c] + M[lr][c-1] + M[lr][c+1]);
				}

				if (R == nrb-2 && rank < nranks-1) {
					isend(&M[rows-2][(C-1)*cbs+1], cbs, rank+1, C, &lowerHaloReqs[C]);
					irecv(&M[rows-1][(C-1)*cbs+1], cbs, rank+1, C, &lowerHaloReqs[C]);
				}
			}
		}
	}
	waitHaloSendRecv(ncb, upperHaloReqs);
	waitHaloSendRecv(ncb, lowerHaloReqs);

	MPI_Barrier(MPI_COMM_WORLD);

	return IGNORE_RESIDUAL;
}
