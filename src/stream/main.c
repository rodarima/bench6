#include "stream.h"

long BS = DEFAULT_BLOCKSIZE;
long SIZE = STREAM_ARRAY_SIZE;

int
main(int argc, char **argv)
{
	(void) argc;
	ssize_t j;
	double t;
	long block;

	if (argc >= 2)
		BS = atol(argv[1]);

	if (argc >= 3)
		SIZE = atol(argv[2]);

	int N_BLOCKS = SIZE/BS;
	if (SIZE%BS != 0) {
		printf("Block size %ld must be divisible by the array size %ld\n",
				BS, SIZE);
		exit(1);
	}

	double *a = malloc(sizeof(double)*(SIZE+OFFSET));
	double *b = malloc(sizeof(double)*(SIZE+OFFSET));
	double *c = malloc(sizeof(double)*(SIZE+OFFSET));

	struct timeval tstart, tend;

	/* Get initial value for system clock. */
	for (block = 0; block < N_BLOCKS; block++) {
		for (j = block*BS; j < (block+1)*BS; j++) {
			a[j] = 1.0;
			b[j] = 2.0;
			c[j] = 0.0;
		}
	}

	t = mysecond();
	for (block = 0; block < N_BLOCKS; block++) {
		for (j = block*BS; j < (block+1)*BS; j++) {
			a[j] = 2.0E0 * a[j];
		}
	}
	t = 1.0E6 * (mysecond() - t);

	gettimeofday(&tstart, NULL);
	kernel_stream(a, b, c, BS, N_BLOCKS);
	gettimeofday(&tend, NULL);
	/*	--- SUMMARY --- */

	/* --- Check Results --- */
	checkSTREAMresults(a, b, c);
	//    printf(HLINE);
	double time_taken = (tend.tv_sec - tstart.tv_sec)*1e6;
	time_taken = (time_taken + (tend.tv_usec - tstart.tv_usec))*1e-6;

	//    printf("Execution time: %lf\n", time_taken);
	printf("%14e %ld %ld %d\n", time_taken, BS, SIZE, NTIMES);

	return 0;
}

/* A gettimeofday routine to give access to the wall
   clock timer on most UNIX-like systems.  */

#include <sys/time.h>

double mysecond(void)
{
        struct timeval tp;

        gettimeofday(&tp,NULL);
        return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6);
}

#ifndef abs
#define abs(a) ((a) >= 0 ? (a) : -(a))
#endif

void
checkSTREAMresults(double *a, double *b, double *c)
{
	STREAM_TYPE aj,bj,cj,scalar;
	STREAM_TYPE aSumErr,bSumErr,cSumErr;
	STREAM_TYPE aAvgErr,bAvgErr,cAvgErr;
	double epsilon;
	ssize_t j;
	int k,ierr;

	/* reproduce initialization */
	aj = 1.0;
	bj = 2.0;
	cj = 0.0;
	/* a[] is modified during timing check */
	aj = 2.0E0 * aj;
	/* now execute timing loop */
	scalar = 3.0;
	for (k=0; k<NTIMES; k++)
	{
		cj = aj;
		bj = scalar*cj;
		cj = aj+bj;
		aj = bj+scalar*cj;
	}

	/* accumulate deltas between observed and expected results */
	aSumErr = 0.0;
	bSumErr = 0.0;
	cSumErr = 0.0;
	for (j=0; j<SIZE; j++) {
		aSumErr += abs(a[j] - aj);
		bSumErr += abs(b[j] - bj);
		cSumErr += abs(c[j] - cj);
		// if (j == 417) printf("Index 417: c[j]: %f, cj: %f\n",c[j],cj);	// MCCALPIN
	}
	aAvgErr = aSumErr / (STREAM_TYPE) SIZE;
	bAvgErr = bSumErr / (STREAM_TYPE) SIZE;
	cAvgErr = cSumErr / (STREAM_TYPE) SIZE;

	if (sizeof(STREAM_TYPE) == 4) {
		epsilon = 1.e-6;
	}
	else if (sizeof(STREAM_TYPE) == 8) {
		epsilon = 1.e-13;
	}
	else {
		printf("WEIRD: sizeof(STREAM_TYPE) = %lu\n",sizeof(STREAM_TYPE));
		epsilon = 1.e-6;
	}

	if (abs(aAvgErr/aj) > epsilon) {
		printf ("Failed Validation on array a[], AvgRelAbsErr > epsilon (%e)\n",epsilon);
		printf ("     Expected Value: %e, AvgAbsErr: %e, AvgRelAbsErr: %e\n",aj,aAvgErr,abs(aAvgErr)/aj);
		ierr = 0;
		for (j=0; j<SIZE; j++) {
			if (abs(a[j]/aj-1.0) > epsilon) {
				ierr++;
			}
		}
		printf("     For array a[], %d errors were found.\n",ierr);
	}
	if (abs(bAvgErr/bj) > epsilon) {
		printf ("Failed Validation on array b[], AvgRelAbsErr > epsilon (%e)\n",epsilon);
		printf ("     Expected Value: %e, AvgAbsErr: %e, AvgRelAbsErr: %e\n",bj,bAvgErr,abs(bAvgErr)/bj);
		printf ("     AvgRelAbsErr > Epsilon (%e)\n",epsilon);
		ierr = 0;
		for (j=0; j<SIZE; j++) {
			if (abs(b[j]/bj-1.0) > epsilon) {
				ierr++;
			}
		}
		printf("     For array b[], %d errors were found.\n",ierr);
	}
	if (abs(cAvgErr/cj) > epsilon) {
		printf ("Failed Validation on array c[], AvgRelAbsErr > epsilon (%e)\n",epsilon);
		printf ("     Expected Value: %e, AvgAbsErr: %e, AvgRelAbsErr: %e\n",cj,cAvgErr,abs(cAvgErr)/cj);
		printf ("     AvgRelAbsErr > Epsilon (%e)\n",epsilon);
		ierr = 0;
		for (j=0; j<SIZE; j++) {
			if (abs(c[j]/cj-1.0) > epsilon) {
				ierr++;
			}
		}
		printf("     For array c[], %d errors were found.\n",ierr);
	}
}
