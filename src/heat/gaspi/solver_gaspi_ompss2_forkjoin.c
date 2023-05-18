#include <GASPI.h>

#include "macros.h"
#include "utils.h"
#include "common/heat.h"


static inline void sendUpperInnerHalo(int64_t rows, int64_t cols, int cbs)
{
	WAIT_FOR_ENTRIES(0, (cols-2)/cbs*2);

	for (int C = 1; C < cols-1; C += cbs) {
		CHECK(gaspi_write_notify(
				UPPER_HALO, (cols+C)*sizeof(double),
				rank-1, LOWER_HALO, (cols+C)*sizeof(double),
				cbs*sizeof(double), C, 1,
				0, GASPI_BLOCK));
	}
}

static inline void sendLowerInnerHalo(int64_t rows, int64_t cols, int cbs)
{
	for (int C = 1; C < cols-1; C += cbs) {
		CHECK(gaspi_write_notify(
				LOWER_HALO, C*sizeof(double),
				rank+1, UPPER_HALO, C*sizeof(double),
				cbs*sizeof(double), C, 1,
				0, GASPI_BLOCK));
	}

	CHECK(gaspi_wait(0, GASPI_BLOCK));
}

static inline void recvUpperOuterHalo(int64_t rows, int64_t cols, int cbs)
{
	gaspi_notification_id_t id;
	gaspi_notification_t val;

	for (int n = 0; n < (cols-2)/cbs; ++n) {
		CHECK(gaspi_notify_waitsome(UPPER_HALO, 0, cols, &id, GASPI_BLOCK));
		CHECK(gaspi_notify_reset(UPPER_HALO, id, &val));
	}
}

static inline void recvLowerOuterHalo(int64_t rows, int64_t cols, int cbs)
{
	gaspi_notification_id_t id;
	gaspi_notification_t val;

	for (int n = 0; n < (cols-2)/cbs; ++n) {
		CHECK(gaspi_notify_waitsome(LOWER_HALO, 0, cols, &id, GASPI_BLOCK));
		CHECK(gaspi_notify_reset(LOWER_HALO, id, &val));
	}
}

static inline void gaussSeidelSolver(int64_t rows, int64_t cols, int rbs, int cbs, double M[rows][cols])
{
	const int nrb = (rows-2)/rbs+2;
	const int ncb = (cols-2)/cbs+2;
	char reps[nrb][ncb];

	if (rank > 0) {
		sendUpperInnerHalo(rows, cols, cbs);
		recvUpperOuterHalo(rows, cols, cbs);
	}
	if (rank < nranks-1)
		recvLowerOuterHalo(rows, cols, cbs);

	for (int R = 1; R < nrb-1; ++R) {
		for (int C = 1; C < ncb-1; ++C) {
			#pragma oss task label("block computation") \
					in(reps[R-1][C]) in(reps[R+1][C]) \
					in(reps[R][C-1]) in(reps[R][C+1]) \
					inout(reps[R][C])
			computeBlock(rows, cols, (R-1)*rbs+1, R*rbs, (C-1)*cbs+1, C*cbs, M);
		}
	}
	#pragma oss taskwait

	if (rank < nranks-1)
		sendLowerInnerHalo(rows, cols, cbs);
}

double solve(HeatConfiguration *conf, int64_t rows, int64_t cols, int timesteps, void *extraData)
{
	double (*matrix)[cols] = (double (*)[cols]) conf->matrix;

	for (int t = 0; t < timesteps; ++t) {
		gaussSeidelSolver(rows, cols, conf->rbs, conf->cbs, matrix);
	}

	CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

	return IGNORE_RESIDUAL;
}
