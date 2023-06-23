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

static char *progname = "b6_runner";

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
	int ret = 0;
	FILE *p = popen(argv[0], "r");

	if (p == NULL) {
		err("popen failed:");
		return -1;
	}

	char line[4096];
	if (fgets(line, 4096, p) == NULL) {
		err("missing stdout line");
		ret = -1;
		goto bad_close;
	}

	char *nl = strchr(line, '\n');
	if (nl != NULL)
		*nl = '\0';

	/* Clean status line */
	fprintf(stderr, "%s\n", line);

	double time;
	sscanf(line, "%le", &time);
	//printf("got %e\n", time);
	*ptime = time;

	/* Drain the rest of the stdout */
	while (fgets(line, 4096, p) != NULL) {
		fprintf(stderr, "%s", line);
	}

bad_close:
	pclose(p);

	return ret;
}

static int
cmp_double(const void *pa, const void *pb)
{
	double a = *(const double *) pa;
	double b = *(const double *) pb;

	if (a < b)
		return -1;
	else if (a > b)
		return 1;
	else
		return 0;
}

static void
stats(struct sampling *s)
{
	if (s->n < 2)
		return;

	double n = s->n;
	//double last = s->samples[s->n - 1];

	/* Sort samples to take the median */
	qsort(s->samples, s->n, sizeof(double), cmp_double);

	double median = s->samples[s->n / 2];

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
	double rstdev = 100.0 * stdev / mean;
	double se = stdev / sqrt(n);
	double rse = 100.0 * se * 1.96 / mean;

	fprintf(stderr, "%s: n=%03d  median=%.3e  mean=%.3e  SD=%.3e  RSD=%.2f%%  RSE=%.2f%%\n",
			progname, s->n, median, mean, stdev, rstdev, rse);

	s->rse = rse;
}

static int
should_continue(struct sampling *s)
{
	stats(s);

	if (s->n < s->nmin)
		return 1;

	if (s->rse > 1.0 /* % */)
		return 1;

	return 0;
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

static int
sample(char *argv[])
{
	struct sampling s = { 0 };
	s.nmax = 4000;
	s.nmin = 30;
	s.samples = calloc(s.nmax, sizeof(double));
	s.n = 0;

	while (should_continue(&s)) {
		double time;
		if (do_run(argv, &time) != 0) {
			err("failed to run benchmark");
			return 1;
		}

		add_sample(&s, time);
	}

	fprintf(stderr, "\n");

	free(s.samples);

	return 0;
}

int
main(int argc, char *argv[])
{
	progname_set(progname);
	(void) argc;

	if (sample(argv+1) != 0) {
		err("failed to sample the benchmark");
		return 1;
	}

	return 0;
}
