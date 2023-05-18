#include <string.h>

#include "common/heat.h"


static inline void gaussSeidelSolver(int64_t rows, int64_t cols, int rbs, int cbs, int nrb, int ncb, double M[rows][cols], char reps[nrb][ncb], double *residual)
{
	for (int R = 1; R < nrb-1; ++R) {
		for (int C = 1; C < ncb-1; ++C) {
			#pragma oss task label("block computation") \
					in(reps[R-1][C]) in(reps[R+1][C]) \
					in(reps[R][C-1]) in(reps[R][C+1]) \
					inout(reps[R][C]) reduction(+: [1]residual)
			*residual += computeBlockResidual(rows, cols, (R-1)*rbs+1, R*rbs, (C-1)*cbs+1, C*cbs, M);
		}
	}
}

double solve(HeatConfiguration *conf, int64_t rows, int64_t cols, int timesteps, void *extraData)
{
	(void) extraData;
	double (*matrix)[cols] = (double (*)[cols]) conf->matrix;
	const double delta = conf->delta;
	const int rbs = conf->rbs;
	const int cbs = conf->cbs;

	const int N = 10;
	double results[N];
	for (int i = 0; i < N; ++i)
		results[i] = delta;

	const int nrb = (rows-2)/rbs+2;
	const int ncb = (cols-2)/cbs+2;
	char representatives[nrb][ncb];

	int t = 0;
	while (t < timesteps) {
		results[t%N] = 0.0f;

		gaussSeidelSolver(rows, cols, rbs, cbs, nrb, ncb, matrix, representatives, &results[t%N]);

		// Advance to the next timestep
		++t;

		// Wait for the Nth previous timestep
		#pragma oss taskwait in(results[t%N])

		// Stop if we reached the residual threshold
		if (results[t%N] < delta)
			break;
	}
	#pragma oss taskwait

	// Save the number of performed timesteps
	conf->convergenceTimesteps = t;

	return results[(t-1)%N];
}
