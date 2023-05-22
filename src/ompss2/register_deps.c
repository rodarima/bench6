/* Copyright (c) 2022 Barcelona Supercomputing Center (BSC)
 * SPDX-License-Identifier: GPL-3.0-or-later */

#include "bench6.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define DEPMAX 500
#define DEPSTEP 10

static long nruns = 30;
static long ntasks = 100;

static int
usage(char *argv[])
{
	fprintf(stderr, "Bench6: A set of Nanos6 micro-benchmarks\n");
	fprintf(stderr, "Usage: %s [-r NRUNS] [-t NTASKS]\n", argv[0]);
	fprintf(stderr, "\n");
	fprintf(stderr,
"Measure the time it takes to create and register NTASKS tasks\n"
"with varying number of dependencies. The number of dependencies\n"
"varies from 0 to 500 in increments of 10. The test is repeated\n"
"NRUNS times.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Defaults: NRUNS=%ld NTASKS=%ld\n", nruns, ntasks);

	return -1;
}

static void
do_run(int run)
{
	int buf[DEPMAX] = { 0 };

	for (int d = 0; d <= DEPMAX; d += DEPSTEP) {
		double t0 = bench6_time();
		for (int t = 0; t < ntasks; t++) {
			#pragma oss task inout({buf[i], i=0;d+1})
			{}
		}
		double t1 = bench6_time();
		printf("%d,%d,%ld,%e\n",
				run, d, ntasks, (t1 - t0) / ntasks);
	}
}

int
main(int argc, char *argv[])
{
	int opt;

	while ((opt = getopt(argc, argv, "hr:t:")) != -1) {
		switch (opt) {
		case 'r':
			nruns = atol(optarg);
			break;
		case 't':
			ntasks = atol(optarg);
			break;
		case 'h': /* Fall through */
		default: /* '?' */
			return usage(argv);
		}
	}

	printf("%s,%s,%s,%s\n", "run", "ndeps", "ntasks", "time_per_task");
	for (int run = 0; run < nruns; run++)
		do_run(run);

	return 0;
}
