#include <TAGASPI.h>

#include "utils.h"
#include "common/streaming.h"

#include <nanos6/debug.h>

static inline void send(int block, int bs, int cbs, int dst, int nqueues)
{
	int queue = nanos6_get_current_virtual_cpu() % nqueues;
	CHECK(tagaspi_write_notify(
			SENDSEG, block*bs*sizeof(double),
			dst, RECVSEG, block*bs*sizeof(double),
			cbs*sizeof(double), block, 1, queue));
}

static inline void recv(int block)
{
	CHECK(tagaspi_notify_async_wait(RECVSEG, block, GASPI_NOTIFICATION_IGNORE));
}

static inline void sendAck(int block, int dst, int nqueues)
{
	int queue = nanos6_get_current_virtual_cpu() % nqueues;
	CHECK(gaspi_notify(SENDSEG, dst, block, 1, queue, GASPI_BLOCK));
}

static inline void ackRemote(int block)
{
	CHECK(tagaspi_notify_async_wait(SENDSEG, block, GASPI_NOTIFICATION_IGNORE));
}

void solve(StreamingConfiguration *conf, uint64_t size, uint64_t bs, int timesteps, void *extraData)
{
	double *array = conf->array;
	double *source = conf->source;
	const uint64_t nb = size/bs + (size%bs > 0);

	const int src = rank - conf->offset;
	const int dst = rank + conf->offset;
	const double compfactor = conf->compfactor;

	GASPInfo *info = (GASPInfo *) extraData;
	const int nqueues = info->nqueues;

	gaspi_notification_t acks[nb];

	for (int t = 0; t < timesteps; ++t) {
		for (uint64_t b = 0; b < nb; ++b) {
			uint64_t cbs = MIN(bs, size-b*bs);
			if (src >= 0) {
				#pragma oss task label("recv") out(source[b*bs;cbs])
				recv(b);
			}

			#pragma oss task label("compute") in(source[b*bs;cbs]) out(array[b*bs;cbs])
			{
				computeBlock(cbs, &array[b*bs], &source[b*bs], compfactor);

				if (t < timesteps-1 && src >= 0) {
					sendAck(b, src, nqueues);
				}
			}

			if (dst < nranks) {
				if (t > 0) {
					#pragma oss task label("ack") out(acks[b])
					ackRemote(b);
				}

				#pragma oss task label("send") in(array[b*bs;cbs]) in(acks[b])
				send(b, bs, cbs, dst, nqueues);
			}
		}
	}
	#pragma oss taskwait

	CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
}
