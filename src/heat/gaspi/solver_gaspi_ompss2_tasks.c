#include <GASPI.h>

#include "macros.h"
#include "utils.h"
#include "common/heat.h"

static int serial;


static inline void sendUpperInnerHalo(int64_t cols, int cbs, int nrb, int ncb, char reps[nrb][ncb])
{
	int coffset = 1;
	for (int C = 1; C < ncb-1; ++C) {
		#pragma oss task label("send upper inner") in(reps[1][C]) inout(serial)
		{
			CHECK(gaspi_write_notify(
					UPPER_HALO, (cols+coffset)*sizeof(double),
					rank-1, LOWER_HALO, (cols+coffset)*sizeof(double),
					cbs*sizeof(double), C, 1,
					0, GASPI_BLOCK));
			CHECK(gaspi_wait(0, GASPI_BLOCK));
		}
		coffset += cbs;
	}
}

static inline void sendLowerInnerHalo(int64_t cols, int cbs, int nrb, int ncb, char reps[nrb][ncb])
{
	int coffset = 1;
	for (int C = 1; C < ncb-1; ++C) {
		#pragma oss task label("send lower inner") in(reps[nrb-2][C]) inout(serial)
		{
			CHECK(gaspi_write_notify(
					LOWER_HALO, coffset*sizeof(double),
					rank+1, UPPER_HALO, coffset*sizeof(double),
					cbs*sizeof(double), C, 1,
					0, GASPI_BLOCK));
			CHECK(gaspi_wait(0, GASPI_BLOCK));
		}
		coffset += cbs;
	}
}

static inline void recvUpperOuterHalo(int nrb, int ncb, char reps[nrb][ncb])
{
	for (int C = 1; C < ncb-1; ++C) {
		#pragma oss task label("recv upper outer") out(reps[0][C]) inout(serial)
		{
			gaspi_notification_id_t id;
			gaspi_notification_t val;
			CHECK(gaspi_notify_waitsome(UPPER_HALO, C, 1, &id, GASPI_BLOCK));
			CHECK(gaspi_notify_reset(UPPER_HALO, C, &val));
		}
	}
}

static inline void recvLowerOuterHalo(int nrb, int ncb, char reps[nrb][ncb])
{
	for (int C = 1; C < ncb-1; ++C) {
		#pragma oss task label("recv lower outer") out(reps[nrb-1][C]) inout(serial)
		{
			gaspi_notification_id_t id;
			gaspi_notification_t val;
			CHECK(gaspi_notify_waitsome(LOWER_HALO, C, 1, &id, GASPI_BLOCK));
			CHECK(gaspi_notify_reset(LOWER_HALO, C, &val));
		}
	}
}

static inline void gaussSeidelSolver(int64_t rows, int64_t cols, int rbs, int cbs, int nrb, int ncb, double M[rows][cols], char reps[nrb][ncb])
{
	if (rank > 0) {
		sendUpperInnerHalo(cols, cbs, nrb, ncb, reps);
		recvUpperOuterHalo(nrb, ncb, reps);
	}
	if (rank < nranks-1)
		recvLowerOuterHalo(nrb, ncb, reps);

	for (int R = 1; R < nrb-1; ++R) {
		for (int C = 1; C < ncb-1; ++C) {
			#pragma oss task label("block computation") \
					in(reps[R-1][C]) in(reps[R+1][C]) \
					in(reps[R][C-1]) in(reps[R][C+1]) \
					inout(reps[R][C])
			computeBlock(rows, cols, (R-1)*rbs+1, R*rbs, (C-1)*cbs+1, C*cbs, M);
		}
	}

	if (rank < nranks-1)
		sendLowerInnerHalo(cols, cbs, nrb, ncb, reps);
}

double solve(HeatConfiguration *conf, int64_t rows, int64_t cols, int timesteps, void *extraData)
{
	double (*matrix)[cols] = (double (*)[cols]) conf->matrix;
	const int rbs = conf->rbs;
	const int cbs = conf->cbs;

	const int nrb = (rows-2)/rbs+2;
	const int ncb = (cols-2)/cbs+2;
	char representatives[nrb][ncb];

	for (int t = 0; t < timesteps; ++t) {
		gaussSeidelSolver(rows, cols, rbs, cbs, nrb, ncb, matrix, representatives);
	}
	#pragma oss taskwait

	CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

	return IGNORE_RESIDUAL;
}
