#include "stream.h"
#include <omp.h>

static void
kernel(double *a, double *b, double *c, int BS, int N_BLOCKS)
{
	int k, j, block;
	double scalar = 3.0;
	for (k = 0; k < NTIMES; k++){
		for (block = 0; block < N_BLOCKS; block++) {
			#pragma omp task \
			depend(out:c[block*BS: (block+1)*BS])\
			depend(in:a[block*BS:(block+1)* BS])
			for (j=block*BS; j<(block+1)*BS; j++)
				c[j] = a[j];
		}
		for (block = 0; block < N_BLOCKS; block++) {
			#pragma omp task \
			depend(out:b[block*BS:(block+1)* BS])\
			depend(in:c[block*BS:(block+1)* BS])	
			for (j=block*BS; j<(block+1)*BS; j++)
				b[j] = scalar*c[j];
		}
		for (block = 0; block < N_BLOCKS; block++) {
			#pragma omp task \
			depend(out:c[block*BS:(block+1)* BS]) \
			depend(in:a[block*BS:(block+1)* BS]) \
			depend(in:b[block*BS:(block+1)* BS])
			for (j=block*BS; j<(block+1)*BS; j++)
				c[j] = a[j]+b[j];
		}
		for (block = 0; block < N_BLOCKS; block++) {
			#pragma omp task \
			depend(out:a[block*BS:(block+1)* BS])\
			depend(in:c[block*BS:(block+1)* BS])\
			depend(in:b[block*BS:(block+1)* BS])
			for (j=block*BS; j<(block+1)*BS; j++)
				a[j] = b[j]+scalar*c[j];
		}
	}
}

void
kernel_stream(double *a, double *b, double *c, int BS, int N_BLOCKS)
{
	#pragma omp parallel
	{
		#pragma omp single
		kernel(a, b, c, BS, N_BLOCKS);
		#pragma omp taskwait
	}
}
