#include <stdint.h>
#include "common/heat.h"

#ifndef SIMD
void computeBlock(const int64_t rows, const int64_t cols,
		const int rstart, const int rend,
		const int cstart, const int cend,
		double M[rows][cols])
{
	for (int r = rstart; r <= rend; ++r) {
		for (int c = cstart; c <= cend; ++c) {
			M[r][c] = 0.25*(M[r-1][c] + M[r+1][c] + M[r][c-1] + M[r][c+1]);
		}
	}
}
#else
void computeBlock(const int64_t rows, const int64_t cols,
		const int rstart, const int rend,
		const int cstart, const int cend,
		double M[rows][cols])
{
	// Assuming square blocks
	const int bs = rend-rstart+1;
	for (int k = 0; k < bs; ++k) {
		#pragma omp simd
		for (int j = 0; j <= k; ++j) {
			const int rr = rstart+k-j;
			const int cc = cstart+j;
			M[rr][cc] = 0.25*(M[rr-1][cc] + M[rr+1][cc] + M[rr][cc-1] + M[rr][cc+1]);
		}
	}
	for (int k = bs-2; k >= 0; --k) {
		#pragma omp simd
		for (int j = 0; j <= k; ++j) {
			const int rr = rstart+bs-j-1;
			const int cc = cstart+bs+j-k-1;
			M[rr][cc] = 0.25*(M[rr-1][cc] + M[rr+1][cc] + M[rr][cc-1] + M[rr][cc+1]);
		}
	}
}
#endif

double computeBlockResidual(const int64_t rows, const int64_t cols,
		const int rstart, const int rend,
		const int cstart, const int cend,
		double M[rows][cols])
{
	double sum = 0.0;
	for (int r = rstart; r <= rend; ++r) {
		for (int c = cstart; c <= cend; ++c) {
			const double value = 0.25*(M[r-1][c] + M[r+1][c] + M[r][c-1] + M[r][c+1]);
			const double diff = value - M[r][c];
			sum += diff*diff;
			M[r][c] = value;
		}
	}
	return sum;
}
