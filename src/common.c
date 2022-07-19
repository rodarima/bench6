/* Copyright (c) 2022 Barcelona Supercomputing Center (BSC)
 * SPDX-License-Identifier: GPL-3.0-or-later */

#define _POSIX_C_SOURCE 199309L

#include "bench6.h"

#include <nanos6/debug.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

/* Returns the current time in seconds since some point in the past */
double get_time()
{
	struct timespec tv;
	if(clock_gettime(CLOCK_MONOTONIC, &tv) != 0)
	{
		perror("clock_gettime failed");
		exit(EXIT_FAILURE);
	}

	return (double)(tv.tv_sec) +
		(double)tv.tv_nsec * 1.0e-9;
}

int get_ncpus()
{
	return (int) nanos6_get_num_cpus();
}
