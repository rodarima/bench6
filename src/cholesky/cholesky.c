#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <cblas.h>
#include <lapacke.h>

#define PRIME1 293
#define PRIME2 719

static void cholesky(long N, long TS, double (*A)[N/TS][TS][TS])
{
	for (long k = 0; k < N/TS; k++) {
		#pragma oss task inout(A[k][k]) label("potrf")
		LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', TS, (double *) A[k][k], TS);
		
		for (long i = k+1; i < N/TS; i++) {
			#pragma oss task in(A[k][k]) inout(A[i][k]) label("trsm")
			cblas_dtrsm(CblasRowMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, TS, TS, 1.0, (double const *) A[k][k], TS, (double *) A[i][k], TS);
		}
		
		for (long i = k+1; i < N/TS; i++) {
			for (long j = k+1; j < i; j++) {
				#pragma oss task in(A[i][k], A[j][k]) inout(A[i][j]) label("gemm")
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, TS, TS, TS, -1.0, (double const *) A[i][k], TS, (double const *) A[j][k], TS, 1.0, (double *) A[i][j], TS);
			}
			#pragma oss task in(A[i][k]) inout(A[i][i]) label("syrk")
			cblas_dsyrk(CblasRowMajor, CblasLower, CblasNoTrans, TS, TS, -1.0, (double const *) A[i][k], TS, 1.0, (double *) A[i][i], TS);
		}
	}
	
	#pragma oss taskwait
}

static void initialize(long N, long TS, double (*a)[N/TS][TS][TS])
{
	for (long i=0; i < N; i++) {
		for (long j=0; j < N; j++) {
			// Generate a value that makes the matrix symmetric positive definite
			double value = (((i+j) % PRIME1) + 1) * (((i+j) % PRIME2) + 1);
			if (i == j) {
				value += PRIME1 * PRIME2;
			}

			// Tiled layout
			a[i/TS][j/TS][i%TS][j%TS] = value;
		}
	}
}

int main(int argc, char **argv)
{
	long n = 4L * 1024L;
	long ts = 1024L;

	if (argc > 3 || (argc >= 2 && strcmp(argv[1], "-h") == 0)) {
		fprintf(stderr, "[USAGE] %s N tasksize\n", argv[0]);
		fprintf(stderr, "  tasksize must divide N\n");
		return 1;
	}

	if (argc >= 2)
		n = atol(argv[1]);

	if (argc >= 3)
		ts = atol(argv[2]);

	if (n % ts != 0) {
		fprintf(stderr, "[USAGE] %s [-h] N tasksize\n", argv[0]);
		fprintf(stderr, "  tasksize must divide N\n");
		return 1;
	}

	typedef double (*matrix_t)[n/ts][ts][ts];

	matrix_t a = (matrix_t) malloc(sizeof(double) * n * n);

	// Warmup iteration
	initialize(n, ts, a);
	cholesky(n, ts, a);

	// Real initialization
	initialize(n, ts, a);

	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC, &start);

	cholesky(n, ts, a);

	clock_gettime(CLOCK_MONOTONIC, &end);

	double duration = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
	duration /= 1000000;

	// GFlops
	double performance = 0.33*n*n*n + 0.5*n*n + 0.17*n;
	performance = performance / duration;
	performance /= 1000000000;

	printf("%14e %14e %14ld %14ld\n", duration, performance, n, ts);

	return 0;
}
