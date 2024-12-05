#include <mpi.h>

#include "utils.h"
#include "common/streaming.h"

typedef struct {
	MPI_Request send;
	MPI_Request recv;
} Requests;

static inline void isend(const double *data, int nelems, int dst, int tag, Requests *req)
{
	MPI_Isend(data, nelems, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &req->send);
}

static inline void irecv(double *data, int nelems, int src, int tag, Requests *req)
{
	MPI_Irecv(data, nelems, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &req->recv);
}

static inline void waitAll(int nreqs, Requests *reqs)
{
	MPI_Waitall(2*nreqs, (MPI_Request *) reqs, MPI_STATUSES_IGNORE);
}

#if 0
static inline void waitSend(Requests *req)
{
	MPI_Wait(&req->send, MPI_STATUS_IGNORE);
}

static inline void waitRecv(Requests *req)
{
	MPI_Wait(&req->recv, MPI_STATUS_IGNORE);
}
#endif

static inline void waitSendRecv(Requests *reqs)
{
	MPI_Waitall(2, (MPI_Request *) reqs, MPI_STATUSES_IGNORE);
}

void solve(StreamingConfiguration *conf, uint64_t size, uint64_t bs, int timesteps, void *extraData)
{
	(void) extraData;
	double *array = conf->array;
	double *source = conf->source;
	const uint64_t nb = size/bs + (size%bs > 0);

	const int src = rank - conf->offset;
	const int dst = rank + conf->offset;
	const double compfactor = conf->compfactor;

	Requests reqs[nb];
	for (uint64_t b = 0; b < nb; ++b) {
		reqs[b].send = MPI_REQUEST_NULL;
		reqs[b].recv = MPI_REQUEST_NULL;
	}

	for (int t = 0; t < timesteps; ++t) {
		if (src >= 0) {
			for (uint64_t b = 0; b < nb; ++b) {
				uint64_t cbs = MIN(bs, size-b*bs);
				irecv(&source[b*bs], cbs, src, b, &reqs[b]);
			}
		}

		for (uint64_t b = 0; b < nb; ++b) {
			uint64_t cbs = MIN(bs, size-b*bs);

			waitSendRecv(&reqs[b]);
			computeBlock(cbs, &array[b*bs], &source[b*bs], compfactor);

			if (dst < nranks) {
				isend(&array[b*bs], cbs, dst, b, &reqs[b]);
			}
		}
	}
	waitAll(nb, reqs);

	MPI_Barrier(MPI_COMM_WORLD);
}
