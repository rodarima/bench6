#include <string.h>
#include <stdio.h>

#include "common/heat.h"

const char *
summary(void)
{
	return "Parallel version using OmpSs-2 tasks and taking into account the\n"
		"residual";
}

static inline void gaussSeidelSolver(int64_t rows, int64_t cols, int rbs, int cbs, int nrb, int ncb, double M[rows][cols], char reps[nrb][ncb], double *residual, double *max_elem, double relax)
{
	for (int R = 1; R < nrb-1; ++R) {
		for (int C = 1; C < ncb-1; ++C) {
			#pragma oss task label("block computation") \
					in(reps[R-1][C]) in(reps[R+1][C]) \
					in(reps[R][C-1]) in(reps[R][C+1]) \
					inout(reps[R][C]) \
					reduction(max: [1]residual) \
					reduction(max: [1]max_elem)
			{
				double lresidual = 0.0;
				double lmax_elem = 0.0;

				computeBlockResidual(rows, cols, (R-1)*rbs+1,
						R*rbs, (C-1)*cbs+1, C*cbs, M,
						relax, &lresidual, &lmax_elem);

				*residual = fmax(*residual, lresidual);
				*max_elem = fmax(*max_elem, lmax_elem);
			}
		}
	}
}

double solve(HeatConfiguration *conf, int64_t rows, int64_t cols, int timesteps, void *extraData)
{
	FILE *f = fopen("convergence.csv", "w");
	fprintf(f, "iter error time\n");

	(void) extraData;
	double (*matrix)[cols] = (double (*)[cols]) conf->matrix;
	const double delta = conf->delta;
	const int rbs = conf->rbs;
	const int cbs = conf->cbs;

	const int N = 4;
	double results[N];
	double max_elem[N];
	double residual[N];

	for (int i = 0; i < N; ++i) {
		results[i] = 666;
		max_elem[i] = 666;
		residual[i] = 666;
	}

	const int nrb = (rows-2)/rbs+2;
	const int ncb = (cols-2)/cbs+2;
	char representatives[nrb][ncb];

	double t0 = getTime();

	int t = 0;
	while (t < timesteps) {
		results[t%N] = 0.0f;
		max_elem[t%N] = 0.0f;
		residual[t%N] = 0.0f;

		gaussSeidelSolver(rows, cols, rbs, cbs, nrb, ncb, matrix,
				representatives, &residual[t%N], &max_elem[t%N], conf->relax);

		#pragma oss task in(residual[t%N], max_elem[t%N]) out(results[t%N])
		{
			results[t%N] = residual[t%N] / max_elem[t%N];
			fprintf(f, "%d %e %e\n", t, results[t%N], getTime() - t0);
			//fprintf(stderr, "t=%d error=%e\n", t, results[t%N]);
		}

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

	fclose(f);

	return results[(t-1)%N];
}
