#ifndef MATMUL_H
#define MATMUL_H

#include <stddef.h>

extern int rank;
extern int nranks;

typedef struct {
	size_t n;
	size_t m;
	size_t ts;
	size_t timesteps;
	size_t warmup;
} matmul_conf_t;

typedef struct {
	void *A;
	void *remote1;
	void *remote2;
	void *B;
	void *C;

	double alpha;
	double beta;
} matmul_t;

double matmul_gettime(void);
void matmul_fail(const char *msg);
void matmul_getconf(int argc, char **argv, matmul_conf_t *conf);
void matmul_setup(matmul_conf_t *conf, matmul_t *matmul, size_t m_per_rank);
void matmul_fill(size_t N, size_t M, size_t TS, double (*A)[M/TS][TS][TS], double value);
int matmul_check(size_t N, size_t M, size_t TS, double (*A)[M/TS][TS][TS], double expected);
void matmul_solve(size_t N, size_t M, size_t TS, matmul_t *mm, size_t timesteps);
void matmul_barrier_issue(size_t N, size_t M, size_t TS, matmul_t *matmul);
void matmul_barrier_notify(void);

#endif
