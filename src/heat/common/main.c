#include <assert.h>
#include <stdio.h>
#include <unistd.h>

#include "common/heat.h"


int main(int argc, char **argv)
{
	HeatConfiguration conf;
	readConfiguration(argc, argv, &conf);
	refineConfiguration(&conf, conf.rbs, conf.cbs);
	if (conf.verbose)
		printConfiguration(&conf);

	int64_t rows = conf.rows+2;
	int64_t cols = conf.cols+2;

	initialize(&conf, rows, cols, 0);

	if (conf.warmup)
		solve(&conf, rows, cols, 1, NULL);

	// Solve the problem
	double start = getTime();
	double residual = solve(&conf, rows, cols, conf.timesteps, NULL);
	double end = getTime();

	int64_t totalElements = conf.rows*conf.cols;
	double throughput = (totalElements*conf.timesteps)/(end-start);

#ifdef _OMPSS_2
	int threads = sysconf(_SC_NPROCESSORS_ONLN);
#else
	int threads = 1;
#endif

	fprintf(stderr,"%8s %8s %8s %8s %8s %8s %14s %14s %14s",
			"rows", "cols", "rbs", "cbs", "threads",
			"steps", "error", "time", "updates/s\n");
	fprintf(stdout, "%8ld %8ld %8d %8d %8d %8d %14e %14e %14e\n", 
			conf.rows, conf.cols, 
			conf.rbs, conf.cbs, threads,
			conf.convergenceTimesteps, residual, end-start, throughput);

	if (conf.generateImage)
		writeImage(conf.imageFileName, conf.matrix, rows, cols);

	finalize(&conf);

	return 0;
}
