#include <mpi.h>
#include <TAMPI.h>

#include "utils.h"
#include "common/streaming.h"

static inline void send(const double *data, int nelems, int dst, int tag)
{
	MPI_Send(data, nelems, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
}

static inline void recv(double *data, int nelems, int src, int tag)
{
	MPI_Recv(data, nelems, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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

	for (int t = 0; t < timesteps; ++t) {
		for (uint64_t b = 0; b < nb; ++b) {
			uint64_t cbs = MIN(bs, size-b*bs);
			if (src >= 0) {
				#pragma oss task label("recv") out(source[b*bs;cbs])
				recv(&source[b*bs], cbs, src, b);
			}

			#pragma oss task label("compute") in(source[b*bs;cbs]) out(array[b*bs;cbs])
			computeBlock(cbs, &array[b*bs], &source[b*bs], compfactor);

			if (dst < nranks) {
				#pragma oss task label("send") in(array[b*bs;cbs])
				send(&array[b*bs], cbs, dst, b);
			}
		}
	}
	#pragma oss taskwait

	MPI_Barrier(MPI_COMM_WORLD);
}
