#include "heat.h"

const char *
summary(void)
{
	return "Sequential solver with one CPU";
}

static inline void gaussSeidelSolver(int64_t rows, int64_t cols, int rbs, int cbs, double M[rows][cols])
{
	for (int R = 1; R < rows-1; R += rbs) {
		for (int C = 1; C < cols-1; C += cbs) {
			computeBlock(rows, cols, R, R+rbs-1, C, C+cbs-1, M);
		}
	}
}

double solve(HeatConfiguration *conf, int64_t rows, int64_t cols, int timesteps, void *extraData)
{
	(void) extraData;
	double (*matrix)[cols] = (double (*)[cols]) conf->matrix;

	for (int t = 0; t < timesteps; ++t) {
		gaussSeidelSolver(rows, cols, conf->rbs, conf->cbs, matrix);
	}

	return IGNORE_RESIDUAL;
}
