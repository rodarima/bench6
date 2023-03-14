/* Copyright (c) 2023 Barcelona Supercomputing Center (BSC)
 * SPDX-License-Identifier: GPL-3.0-or-later */

#define _DEFAULT_SOURCE

#include "bench6.h"

#include <nanos6.h>
#include <nanos6/debug.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

static int ncpus = -1;
static long nruns = 300L;
static long ntasks_per_cpu = 5000L;
static double taskwork_us = 10.0;

static double t0;

#define M_WORK 10000000L

static void
busywork(long loops)
{
	for (volatile long j = 0; j < loops; j++);
}

static double
get_time_ms(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double) ts.tv_sec + (double) ts.tv_nsec * 1.0e-9;
}

static void
dummy_work(double ms)
{
	double end = get_time_ms() + ms * 1e-3;
	while (get_time_ms() < end);
}

static void
do_run(int run)
{
	/* Warm up all the threads */
	for (long i = 0L; i < ncpus; i++) {
		#pragma oss task label("warmup")
		dummy_work(20.0);
	}

	#pragma oss taskwait

	int flag = 0;

	/* Delay task start so we can create lots of tasks */
	#pragma oss task inout(flag) label("trigger")
	{
		dummy_work(10.0);
		busywork(M_WORK); /* So we can see it in perf */
		t0 = get_time();
	}

	for (long i = 0L; i < ntasks_per_cpu * ncpus; i++) {
		#pragma oss task in(flag) label("quickie")
		dummy_work(taskwork_us * 1e-3);
	}

	/* When trigger finishes all small tasks will become ready */
	#pragma oss taskwait
	double t1 = get_time();
	printf("%d,%ld,%d,%.3f,%.3f\n", run, ntasks_per_cpu, ncpus, taskwork_us, (t1 - t0) * 1e3);
}

static int
usage(char *argv[])
{
	fprintf(stderr, "Bench6: A set of Nanos6 micro-benchmarks\n");
	fprintf(stderr, "Usage: %s [-r NRUNS] [-t NTASKS_PER_CPU] [-w TASKWORK_US]\n", argv[0]);
	fprintf(stderr, "\n");
	fprintf(stderr,
"Creates NTASKS_PER_CPU tasks per CPU that become ready at the\n"
"same time, like a wave. The time between the moment they become\n"
"ready and when the all end is measured. Tasks run for TASKWORK_US\n"
"microseconds. The test is repeated NRUNS times.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Defaults: NRUNS=%ld NTASKS_PER_CPU=%ld TASKWORK_US=%.3f\n",
			nruns, ntasks_per_cpu, taskwork_us);

	return -1;
}

int
main(int argc, char *argv[])
{
	int opt;

	while ((opt = getopt(argc, argv, "hr:t:w:")) != -1) {
		switch (opt) {
		case 'r':
			nruns = atol(optarg);
			break;
		case 't':
			ntasks_per_cpu = atol(optarg);
			break;
		case 'w':
			taskwork_us = atof(optarg);
			break;
		case 'h': /* Fall through */
		default: /* '?' */
			return usage(argv);
		}
	}

	ncpus = get_ncpus();

	printf("%s,%s,%s,%s,%s\n", "run", "ntasks_per_cpu", "ncpus", "taskwork_us", "time_ms");
	for (int run = 0; run < nruns; run++)
		do_run(run);

	return 0;
}
