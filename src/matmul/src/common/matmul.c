#include "matmul.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double matmul_gettime(void)
{
	struct timespec tv;
	clock_gettime(CLOCK_MONOTONIC, &tv);
	return tv.tv_sec+1e-9*tv.tv_nsec;
}

void matmul_fail(const char *msg)
{
	fprintf(stderr, "error: %s\n", msg);
	exit(1);
}

void matmul_getconf(int argc, char **argv, matmul_conf_t *conf)
{
	if (argc != 5) {
		fprintf(stderr, "Usage: %s N TIMESTEPS TILESIZE WARMUP", argv[0]);
		matmul_fail("missing parameters");
	}

	conf->n = conf->m = atoi(argv[1]);
	conf->timesteps = atoi(argv[2]);
	conf->ts = atoi(argv[3]);
	conf->warmup = atoi(argv[4]);
}

void matmul_setup(matmul_conf_t *conf, matmul_t *matmul, size_t m_per_rank)
{
	const size_t n = conf->n;
	const size_t m = m_per_rank;
	const size_t ts = conf->ts;

	const size_t matrix_size = n * m * sizeof(double);

	typedef double (*matrix_MxN_t)[n/ts][ts][ts];
	typedef double (*matrix_NxM_t)[m/ts][ts][ts];

	matrix_MxN_t A, remote1, remote2;
	matrix_NxM_t B, C;
	A = malloc(matrix_size);
	B = malloc(matrix_size);
	C = malloc(matrix_size);
	remote1 = malloc(matrix_size);
	remote2 = malloc(matrix_size);

	if (!A || !B || !C || !remote1 || !remote2)
		matmul_fail("not enough memory");

	// Initialize matrices
	matmul_fill(m, n, ts, A, 1.0);
	matmul_fill(n, m, ts, B, 1.0);
	matmul_fill(n, m, ts, C, 0.0);

	matmul->A = (void *)A;
	matmul->remote1 = (void *)remote1;
	matmul->remote2 = (void *)remote2;
	matmul->B = (void *)B;
	matmul->C = (void *)C;
	matmul->alpha = 1.0;
	matmul->beta = 1.0;
}

void matmul_fill(size_t N, size_t M, size_t TS, double (*A)[M/TS][TS][TS], double value)
{
	for (size_t R = 0; R < N/TS; ++R) {
		for (size_t C = 0; C < M/TS; ++C) {
			for (size_t r = 0; r < TS; ++r) {
				for (size_t c = 0; c < TS; ++c) A[R][C][r][c] = value;
			}
		}
	}
}

int matmul_check(size_t N, size_t M, size_t TS, double (*A)[M/TS][TS][TS], double expected)
{
	int errors = 0;
	for (size_t R = 0; R < N/TS; ++R) {
		for (size_t C = 0; C < M/TS; ++C) {
			for (size_t r = 0; r < TS; ++r) {
				for (size_t c = 0; c < TS; ++c) {
					double value = A[R][C][r][c];
					if (fabs(value-expected) > 1e-5) {
						++errors;
						fprintf(stderr, "C[%ld][%ld][%ld][%ld]: %f, expected: %f\n", R, C, r, c, value, expected);
					}
				}
			}
		}
	}
	return errors;
}
