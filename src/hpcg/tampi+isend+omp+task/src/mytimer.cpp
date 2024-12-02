
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/////////////////////////////////////////////////////////////////////////

// Function to return time in seconds.
// If compiled with no flags, return CPU time (user and system).
// If compiled with -DWALL, returns elapsed time.

/////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <time.h>
#include <sys/resource.h>

double mytimer(void);

double mytimer(void) {
  struct timespec tp;
  static long start=0, startn;
  if (!start) {
    clock_gettime(CLOCK_MONOTONIC, &tp);
    start = tp.tv_sec;
    startn = tp.tv_nsec;
    return 0.0;
  }
  clock_gettime(CLOCK_MONOTONIC, &tp);
  return ((double) (tp.tv_sec - start)) + (tp.tv_nsec-startn)/1000000000.0 ;
}
