/* Copyright (c) 2022 Barcelona Supercomputing Center (BSC)
 * SPDX-License-Identifier: GPL-3.0-or-later */

#include "bench6.h"

#include <stdio.h>
#include <string.h>

struct bench {
	int (*fn) (int, char *[]);
	char *name;
};

struct bench benchmarks[] = {
	{ bench6_sched_get,     "sched_get" },
	{ bench6_sched_add,     "sched_add" },
	{ bench6_register_deps, "register_deps" },
	{ NULL, NULL }
};

static int
usage()
{
	fprintf(stderr, "Bench6: A set of Nanos6 micro-benchmarks\n");
	fprintf(stderr, "Usage: bench6 BENCHMARK [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Available benchmarks:\n");
	for (struct bench *p = &benchmarks[0]; p->name != NULL; p++)
		fprintf(stderr, "  %s\n", p->name);
	fprintf(stderr, "\n");
	fprintf(stderr, "Use \"bench6 BENCHMARK -h\" for specific options.\n");

	return -1;
}

int
main(int argc, char *argv[])
{
	if (argc <= 1)
		return usage();

	char *name = argv[1];

	argc--;
	argv++;

	for (struct bench *p = &benchmarks[0]; p->name != NULL; p++) {
		if (strcmp(p->name, name) == 0)
			return p->fn(argc, argv);
	}

	return usage();
}
