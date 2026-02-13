#include <stdatomic.h>
#include <assert.h>
#include <math.h>     /* pow */
#include <stdio.h>    /* printf */
#include <stdlib.h>
#include <sys/time.h> /* gettimeofday */
#include <string.h>

#ifdef NO_TASKWAIT
atomic_long results[4];
atomic_long final_result;
#else
long results[4];
long final_result;
#endif

long size;
long final_depth;
size_t NUM_ITER = 4;
double max_seconds;
struct timeval start;

#ifdef NO_TASKWAIT
static void fibonacci(long index, atomic_long *resultPointer)
{
	if (index <= 1) {
		atomic_fetch_add(resultPointer, index);
	} else {
		#pragma oss task label("fibonacci") /* cost(pow(2,index-1)) */ final((index-1) <= final_depth)
		fibonacci(index-1, resultPointer);

		#pragma oss task label("fibonacci") /* cost(pow(2,index-2)) */ final((index-2) <= final_depth)
		fibonacci(index-2, resultPointer);
	}
}
#else
static void fibonacci(long index, long *resultPointer)
{
	if (index <= 1) {
		*resultPointer = index;
		return;
	}

	long result1, result2;

	#pragma oss task label("fibonacci") shared(result1) /* cost(pow(2,index-1)) */ final((index-1) <= final_depth)
	fibonacci(index-1, &result1);

	#pragma oss task label("fibonacci") shared(result2) /* cost(pow(2,index-2)) */ final((index-2) <= final_depth)
	fibonacci(index-2, &result2);

	#pragma oss taskwait
	*resultPointer = result1 + result2;
}
#endif

static void oss_fibonacci(void)
{
	if (max_seconds != 0) {
		size_t it       = 0;
		double elapsed  = 0;
		struct timeval stop;

		gettimeofday(&stop, NULL);
		elapsed = 1000000 * (stop.tv_sec - start.tv_sec);
		elapsed += stop.tv_usec - start.tv_usec;
		elapsed /= 1000000.0;
		while (elapsed < max_seconds) {
			results[it % NUM_ITER] = 0;

			#pragma oss task label("fibonacci") \
				/* cost(pow(2,size)) */       \
				final(size <= final_depth)    \
				out(results[it % NUM_ITER])
			fibonacci(size, &results[it % NUM_ITER]);

			#pragma oss taskwait in(results[(it+1) % NUM_ITER])

			gettimeofday(&stop, NULL);
			elapsed = 1000000 * (stop.tv_sec - start.tv_sec);
			elapsed += stop.tv_usec - start.tv_usec;
			elapsed /= 1000000.0;

			it++;
		}

		final_result = results[(it-1) % NUM_ITER];
	} else {
		#pragma oss task label("fibonacci") shared(final_result) /* cost(pow(2,size)) */ final(size <= final_depth)
		fibonacci(size, &final_result);
	}
}


int main(int argc, char **argv)
{
	if (argc > 4 || (argc >= 2 && strcmp(argv[1], "-h") == 0)) {
		printf("\nUsage: %s [-h] <N> <F> [S]\n", argv[0]);
		printf("  N : The size of the problem / how many numbers must be generated (fibonacci sequence size)\n");
		printf("  F : The maximum depth (numbers generated) in which tasks become final\n");
		printf("  S : Optional, minimum amount of seconds the execution must take (as many repetitions as needed\n\n");
		return 1;
	}

	// Initialize
	size = 20;
	final_depth = 5;
	max_seconds  = 0.0;
	final_result = 0;

	if (argc >= 2)
		size = atoi(argv[1]);
	if (argc >= 3)
		final_depth = atoi(argv[2]);

	if (argc >= 4) {
		max_seconds = atof(argv[3]);
		assert(max_seconds > 0);
	}

	// Start timing
	gettimeofday(&start, NULL);

	// Compute the sequence as many times as needed
	oss_fibonacci();
	#pragma oss taskwait

	// Stop timing
	struct timeval stop;
	gettimeofday(&stop, NULL);
	unsigned long elapsed = 1000000 * (stop.tv_sec - start.tv_sec);
	elapsed += stop.tv_usec - start.tv_usec;

	double t = (double) elapsed / 1000000.0;

	printf("%14e %zu %zu %zu %s\n", t, size, final_depth, final_result, BENCH6_NAME);

	// Print results
	//printf("\n");
	//printf("%s: %s\n" , "BENCHMARK"  , "fibonacci");
	//printf("%s: %zu\n", "RESULT"     , final_result);
	//printf("%s: %zu\n", "SIZE"       , size);
	//printf("%s: %zu\n", "FINAL_DEPTH", final_depth);
	//printf("%s: %f\n" , "TIME(s)"    , (double)elapsed/1000000.0);
	//printf("\n");

}
