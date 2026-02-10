#include "stream.h"

static void
kernel(double *a, double *b, double *c, int BS, int N_BLOCKS)
{
	int k, j, block;
	double scalar = 3.0;
	for (k = 0; k < NTIMES; k++) {
		for (block = 0; block < N_BLOCKS; block++) {
			#pragma oss task out(c[block*BS]) in(a[block*BS])
			for (j=block*BS; j<(block+1)*BS; j++)
				c[j] = a[j];
		}
		for (block = 0; block < N_BLOCKS; block++) {
			#pragma oss task out(b[block*BS]) in(c[block*BS])
			for (j=block*BS; j<(block+1)*BS; j++)
				b[j] = scalar*c[j];
		}
		for (block = 0; block < N_BLOCKS; block++) {
			#pragma oss task out(c[block*BS]) in(a[block*BS]) in(b[block*BS])
			for (j=block*BS; j<(block+1)*BS; j++)
				c[j] = a[j]+b[j];
		}
		for (block = 0; block < N_BLOCKS; block++) {
			#pragma oss task out(a[block*BS]) in(c[block*BS]) in(b[block*BS])
			for (j=block*BS; j<(block+1)*BS; j++)
				a[j] = b[j]+scalar*c[j];
		}
	}
}

void
kernel_stream(double *a, double *b, double *c, int BS, int N_BLOCKS)
{
	kernel(a, b, c, BS, N_BLOCKS);
	#pragma oss taskwait
}
