/* Copyright (c) 2022 Barcelona Supercomputing Center (BSC)
 * SPDX-License-Identifier: GPL-3.0-or-later */

#include "bench6.h"

#include <nanos6.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

static int ncpus = -1;
static long nruns = 30L;
static long ntasks = 10000L;

static atomic_int wait = 0;
static void **handle;

static int
usage(char *argv[])
{
	fprintf(stderr, "Bench6: A set of Nanos6 micro-benchmarks\n");
	fprintf(stderr, "Usage: %s [-r NRUNS] [-t NTASKS]\n", argv[0]);
	fprintf(stderr, "\n");
	fprintf(stderr,
"Measure the time it takes to unblock NTASKS which will end\n"
"immediately. The test is repeated NRUNS times.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Defaults: NRUNS=%ld NTASKS=%ld\n", nruns, ntasks);

	return -1;
}

#pragma oss task
static void
do_task(int t)
{
	if (atomic_load(&wait)) {
		handle[t] = nanos6_get_current_blocking_context();
		nanos6_block_current_task(handle[t]);
	}
}

static void
do_run(int run)
{
	memset(handle, 0, ntasks * sizeof(void *));
	atomic_fetch_add(&wait, 1);

	for (int t = 0; t < ntasks; t++)
		do_task(t);

	double t0 = get_time();
	atomic_fetch_sub(&wait, 1);

	for (int t = 0; t < ntasks; t++) {
		if (handle[t]) {
			nanos6_unblock_task(handle[t]);
		}
	}

	#pragma oss taskwait

	double t1 = get_time();
	printf("%d,%ld,%e,%e\n",
			run, ntasks, (t1 - t0),
			(t1 - t0) / ((double) ntasks) * ((double) ncpus));
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

	ncpus = get_ncpus();

	handle = calloc(ntasks, sizeof(void *));

	if (handle == NULL) {
		perror("calloc failed");
		return -1;
	}

	printf("%s,%s,%s,%s\n", "run", "ntasks", "time", "time_per_task_per_cpu");
	for (int run = 0; run < nruns; run++)
		do_run(run);

	return 0;
}
