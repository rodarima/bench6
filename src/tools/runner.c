#include "common.h"
#include "config.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

//static void
//usage(void)
//{
//	exit(1);
//}

static const char *cleanup_cmd = NULL;
static char **bench_argv = NULL;

struct sampling {
	int nmax;
	int nmin;
	int n;
	double *samples;
	double rse;
	double last;
};

static int
do_run(char *argv[], double *ptime)
{
	/* Gather binary path */
	char path[PATH_MAX];
	sprintf(path, "%s/%s", BENCH6_BIN, argv[0]);

	if (access(path, R_OK | X_OK) != 0) {
		err("cannot find benchmark %s:", path);
		return -1;
	}

	int pipefd[2];
	if (pipe(pipefd) != 0) {
		err("pipe failed:");
		return -1;
	}

	/* Fork */
	pid_t p = fork();

	if (p < 0) {
		err("fork failed:");
		return -1;
	}

	/* In children execute benchmark */
	if (p == 0) {
		close(pipefd[0]);
		dup2(pipefd[1], 1);
		close(2);
		if (execve(path, argv, NULL) != 0) {
			err("execve failed:");
			return -1;
		}
		/* Not reached */
	} else {
		close(pipefd[1]);
		char line[4096];
		FILE *f = fdopen(pipefd[0], "r");
		if (f == NULL) {
			err("fdopen failed:");
			return -1;
		}

		if (fgets(line, 4096, f) == NULL) {
			err("missing stdout line");
			return -1;
		}

		char *nl = strchr(line, '\n');
		if (nl != NULL)
			*nl = '\0';

		double time;
		sscanf(line, "%le", &time);
		//printf("got %e\n", time);
		*ptime = time;

		/* Drain the rest of the stdout */
		while (fgets(line, 4096, f) != NULL) { }
		fclose(f);
		close(pipefd[0]);
	}

	if (cleanup_cmd != NULL)
		system(cleanup_cmd);

	return 0;
}

static void
stats(struct sampling *s)
{
	if (s->n < 2)
		return;

	double n = s->n;
	double sum = 0.0;
	for (int i = 0; i < s->n; i++)
		sum += s->samples[i];

	double mean = sum / n;
	double sumsqr = 0.0;
	for (int i = 0; i < s->n; i++) {
		double dev = s->samples[i] - mean;
		sumsqr += dev * dev;
	}

	double var = sumsqr / n;
	double stdev = sqrt(var);
	double se = stdev / sqrt(n);
	double rse = se * 1.96 / mean;

	fprintf(stderr, "n=%d last=%e mean=%e stdev=%e se=%e rse=%e\n",
			s->n, s->last, mean, stdev, se, rse);

	s->rse = rse;
}

static int
should_stop(struct sampling *s)
{
	stats(s);

	if (s->n < s->nmin)
		return 0;

	if (s->rse * 100.0 < 1.0 /* % */)
		return 0;

	return 1;
}

static void
add_sample(struct sampling *s, double time)
{
	if (s->n >= s->nmax) {
		die("overflowing samples");
	} else {
		s->samples[s->n] = time;
		s->n++;
		s->last = time;
	}
}

//static int
//compare_double(const void *a, const void *b)
//{
//	double aa = *(const double *) a;
//	double bb = *(const double *) b;
//
//	if (aa < bb)
//		return -1;
//	else if (aa > bb)
//		return +1;
//	else
//		return 0;
//}

static int
sample(char *argv[])
{
	struct sampling s = { 0 };
	s.nmax = 4000;
	s.nmin = 30;
	s.samples = calloc(s.nmax, sizeof(double));
	s.n = 0;

	while (!should_stop(&s)) {
		double time;
		if (do_run(argv, &time) != 0) {
			err("failed to run benchmark");
			return 1;
		}

		add_sample(&s, time);
	}

	free(s.samples);

	return 0;
}

static void
usage(void)
{
	fprintf(stderr, "c:h\n");
	exit(1);
}

static void
parse_args(int argc, char *argv[])
{
	int opt;

	while ((opt = getopt(argc, argv, "c:h")) != -1) {
		switch (opt) {
			case 'c':
				cleanup_cmd = optarg;
				break;
			case 'h':
			default: /* '?' */
				usage();
		}
	}

	if (optind >= argc) {
		err("bad usage: program");
		usage();
	}

	bench_argv = &argv[optind];
}

int
main(int argc, char *argv[])
{
	(void) argc;

	parse_args(argc, argv);

	if (sample(bench_argv) != 0) {
		err("failed to sample the benchmark");
		return 1;
	}

	return 0;
}
