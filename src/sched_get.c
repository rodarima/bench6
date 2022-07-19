/* Copyright (c) 2022 Barcelona Supercomputing Center (BSC)
 * SPDX-License-Identifier: GPL-3.0-or-later */

#define _POSIX_C_SOURCE 2

#include "bench6.h"

#include <stdio.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <unistd.h>

static long nruns = 30L;
static long ntasks = 20000L;

static atomic_int wait = 0;

static int
usage(char *argv[])
{
	fprintf(stderr, "Usage: bench6 %s [-r NRUNS] [-t NTASKS]\n", argv[0]);
	fprintf(stderr, "\n");
	fprintf(stderr,
"Creates NTASKS tasks without dependencies, but the tasks don't\n"
"finish until the creator ends. As soon as all tasks are created\n"
"the time is measured until all tasks end. The test is repeated\n"
"NRUNS times.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Defaults: NRUNS=%ld NTASKS=%ld\n", nruns, ntasks);

	return -1;
}

static void
do_run(int run)
{
	atomic_fetch_add(&wait, 1);

	for (int t = 0; t < ntasks; t++) {
		#pragma oss task
		{
			while (atomic_load(&wait));
		}
	}

	double t0 = get_time();
	atomic_fetch_sub(&wait, 1);

	#pragma oss taskwait
	double t1 = get_time();
	printf("%d,%ld,%e,%e\n",
			run, ntasks, (t1 - t0),
			(t1 - t0) / ((double) ntasks));
}

int
bench6_sched_get(int argc, char *argv[])
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

	printf("%s,%s,%s,%s\n", "run", "ntasks", "time", "time_per_task");
	for (int run = 0; run < nruns; run++)
		do_run(run);

	return 0;
}
