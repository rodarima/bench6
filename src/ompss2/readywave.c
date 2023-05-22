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
#include <stdatomic.h>
#include <math.h>
#include <pthread.h>

static char progname[] = "bench6.readywave";
static int ncpus = -1;
static long nwarm = 100L;
static long nruns = 200L;
static long ntasks_per_cpu = 1000L;
static double size_per_cpu_ns = 400.0;
static double cooldown_ms = 0.0;
static int sequential_sched = 0;

static atomic_int wait = 0;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

#define M_WORK 10000000L

static void
busywork(long loops)
{
	for (volatile long j = 0; j < loops; j++);
}

static double
dummy_work(double ms)
{
	double start = bench6_time();
	double end = start + ms * 1e-3;
	double last;
	while ((last = bench6_time()) < end) {
		busywork(100L);
	}

	return (last - start) * 1e3;
}

static void
do_run(int run)
{
	dummy_work(cooldown_ms);

	atomic_store(&wait, 1);

	/* Cover 2 times the number of CPUs so no quickie starts before the
	 * blockers */
	for (long i = 0L; i < 2*ncpus; i++) {
		#pragma oss task label("blocker")
		{
			//fprintf(stderr, "blocker %d up\n", i);
			/* Wait until the creator finishes */
			while (atomic_load(&wait));
		}
	}

	/* Create the quick tasks */
	for (long i = 0L; i < ntasks_per_cpu * ncpus; i++) {
		#pragma oss task label("quickie")
		{
			if (sequential_sched) {
				dummy_work((ncpus - 1) * size_per_cpu_ns * 1e-6);

				pthread_mutex_lock(&mutex);
				dummy_work(size_per_cpu_ns * 1e-6);
				pthread_mutex_unlock(&mutex);
			} else {
				dummy_work(ncpus * size_per_cpu_ns * 1e-6);
			}
		}
	}

	/* Release the blockers */
	atomic_fetch_sub(&wait, 1);

	/* Start counting the time as the quickies will run now */
	double t0 = bench6_time();

	/* Wait until all tasks are ready */
	#pragma oss taskwait

	/* And measure the end time */
	double t1 = bench6_time();

	/* Warmup run */
	if (run < 0)
		return;

	printf("%d,%ld,%d,%.3f,%e,%e\n",
			run, ntasks_per_cpu,
			ncpus, size_per_cpu_ns,
			(t1 - t0) * 1e3,
			(t1 - t0) * 1e9 / ntasks_per_cpu / ncpus);
}

static void
do_warmup(void)
{
	fprintf(stderr, "running %ld warmup iterations...\n", nwarm);
	/* Warm up all the workers */
	for (long i = 0L; i < 5*ncpus; i++) {
		#pragma oss task label("warmup")
		dummy_work(20.0);
	}

	#pragma oss taskwait

	for (int i = 0; i < nwarm; i++)
		do_run(-1);

	#pragma oss taskwait
	fprintf(stderr, "warmup done\n");
}

static int
usage(void)
{
	fprintf(stderr, "%s - Create a wave of ready rasks\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: %s [-w NWARM] [-r NRUNS] [-t NTASKS] [-s SIZE]\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Creates a large number of ready tasks to put pressure in the\n");
	fprintf(stderr, "scheduler server. First, 2*ncpus tasks block the cpus with\n");
	fprintf(stderr, "work until the creator worker finishes creating all the tasks.\n");
	fprintf(stderr, "Then, the blocker tasks are signaled to finish, and the quickie\n");
	fprintf(stderr, "tasks follow. The time is measured from the signal until they\n");
	fprintf(stderr, "all end.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -w   Number of warmup repetitions. These are used to remove\n");
	fprintf(stderr, "       the effect of the jemalloc contention while growing the\n");
	fprintf(stderr, "       arenas (default %ld).\n", nwarm);
	fprintf(stderr, "\n");
	fprintf(stderr, "  -r   Number of repetitions of the test (default %ld).\n", nruns);
	fprintf(stderr, "\n");
	fprintf(stderr, "  -t   Number of tasks per CPU to be created (default %ld).\n", ntasks_per_cpu);
	fprintf(stderr, "\n");
	fprintf(stderr, "  -s   Size of the tasks in ns per CPU (default %f).\n", size_per_cpu_ns);
	fprintf(stderr, "\n");
	fprintf(stderr, "  -c   Cooldown delay in milliseconds before a new run (default %f).\n", cooldown_ms);
	fprintf(stderr, "\n");
	fprintf(stderr, "  -S   Serve the tasks sequentially (default %s).\n",
			sequential_sched ? "yes" : "no");
	fprintf(stderr, "\n");

	return -1;
}

int
main(int argc, char *argv[])
{
	int opt;

	while ((opt = getopt(argc, argv, "hr:w:t:s:c:S")) != -1) {
		switch (opt) {
		case 'r':
			nruns = atol(optarg);
			break;
		case 'w':
			nwarm = atol(optarg);
			break;
		case 't':
			ntasks_per_cpu = atol(optarg);
			break;
		case 's':
			size_per_cpu_ns = atof(optarg);
			break;
		case 'c':
			cooldown_ms = atof(optarg);
			break;
		case 'S':
			sequential_sched = 1;
			break;
		case 'h': /* Fall through */
		default: /* '?' */
			return usage();
		}
	}

	ncpus = get_ncpus();

	printf("%s,%s,%s,%s,%s,%s\n", "run", "ntasks_per_cpu", "ncpus", "size_per_cpu_ns", "time_ms", "avg_serve_time_ns");
	do_warmup();

	for (int run = 0; run < nruns; run++)
		do_run(run);

	return 0;
}
