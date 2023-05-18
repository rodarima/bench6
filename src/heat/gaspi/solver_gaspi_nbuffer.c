#include <GASPI.h>

#include "macros.h"
#include "utils.h"
#include "common/heat.h"


static inline void startUpperInnerHaloSend(int64_t rows, int64_t cols, int cbs, int C, const GASPInfo *info)
{
	const gaspi_queue_id_t queue = info->queues[C];
	const int coffset = cols+(C-1)*cbs+1;

	WAIT_FOR_ENTRIES(queue, 2);

	CHECK(gaspi_write_notify(
			UPPER_HALO, coffset*sizeof(double),
			rank-1, LOWER_HALO, coffset*sizeof(double),
			cbs*sizeof(double), C, 1,
			queue, GASPI_BLOCK));
}

static inline void startLowerInnerHaloSend(int64_t rows, int64_t cols, int cbs, int C, const GASPInfo *info)
{
	const gaspi_queue_id_t queue = info->queues[C];
	const int coffset = (C-1)*cbs+1;

	WAIT_FOR_ENTRIES(queue, 2);

	CHECK(gaspi_write_notify(
			LOWER_HALO, coffset*sizeof(double),
			rank+1, UPPER_HALO, coffset*sizeof(double),
			cbs*sizeof(double), C, 1,
			queue, GASPI_BLOCK));
}

static inline void waitInnerHaloSend(int C, const GASPInfo *info)
{
	CHECK(gaspi_wait(info->queues[C], GASPI_BLOCK));
}

static inline void waitUpperOuterHaloRecv(int C)
{
	gaspi_notification_id_t id;
	gaspi_notification_t val;
	CHECK(gaspi_notify_waitsome(UPPER_HALO, C, 1, &id, GASPI_BLOCK));
	CHECK(gaspi_notify_reset(UPPER_HALO, C, &val));
}

static inline void waitLowerOuterHaloRecv(int C)
{
	gaspi_notification_id_t id;
	gaspi_notification_t val;
	CHECK(gaspi_notify_waitsome(LOWER_HALO, C, 1, &id, GASPI_BLOCK));
	CHECK(gaspi_notify_reset(LOWER_HALO, C, &val));
}

double solve(HeatConfiguration *conf, int64_t rows, int64_t cols, int timesteps, void *extraData)
{
	double (*M)[cols] = (double (*)[cols]) conf->matrix;
	const GASPInfo *info = (GASPInfo *)extraData;
	const int rbs = conf->rbs;
	const int cbs = conf->cbs;
	const int nrb = (rows-2)/rbs+2;
	const int ncb = (cols-2)/cbs+2;

	for (int t = 0; t < timesteps; ++t) {
		for (int R = 1; R < nrb-1; ++R) {
			for (int C = 1; C < ncb-1; ++C) {
				if (R == 1 && rank > 0)
					waitUpperOuterHaloRecv(C);

				// First row of the block
				const int fr = (R-1)*rbs+1;
				for (int c = (C-1)*cbs+1; c <= C*cbs; ++c) {
					M[fr][c] = 0.25*(M[fr-1][c] + M[fr+1][c] + M[fr][c-1] + M[fr][c+1]);
				}

				if (fr == 1 && rank > 0)
					startUpperInnerHaloSend(rows, cols, cbs, C, info);

				computeBlock(rows, cols, (R-1)*rbs+2, R*rbs-1, (C-1)*cbs+1, C*cbs, M);

				const int lr = R*rbs;
				if (lr == rows-2 && rank < nranks-1) {
					waitInnerHaloSend(C, info);
					if (t > 0) waitLowerOuterHaloRecv(C);
				}

				// Last row of the block
				for (int c = (C-1)*cbs+1; c <= C*cbs; ++c) {
					M[lr][c] = 0.25*(M[lr-1][c] + M[lr+1][c] + M[lr][c-1] + M[lr][c+1]);
				}

				if (R == nrb-2 && rank < nranks-1)
					startLowerInnerHaloSend(rows, cols, cbs, C, info);
			}
		}

		if (t == timesteps-1 && rank < nranks-1) {
			for (int C = 1; C < ncb-1; ++C) {
				waitInnerHaloSend(C, info);
				waitLowerOuterHaloRecv(C);
			}
		}
	}

	CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

	return IGNORE_RESIDUAL;
}
