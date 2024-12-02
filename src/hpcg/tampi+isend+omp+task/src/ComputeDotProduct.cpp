
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

/*!
 @file ComputeDotProduct.cpp

 HPCG routine
 */

#include "mytimer.hpp"
#include "ComputeDotProduct.hpp"
#include "ComputeDotProduct_ref.hpp"

/*!
  Routine to compute the dot product of two vectors.

  This routine calls the reference dot-product implementation by default, but
  can be replaced by a custom routine that is optimized and better suited for
  the target system.

  @param[in] [start, end) the range of vector elements (on this processor)
  @param[in]  x, y the input vectors
  @param[out] result a pointer to scalar value, on exit will contain the result.
  @param[out] time_allreduce the time it took to perform the communication between processes
  @param[out] isOptimized should be set to false if this routine uses the reference implementation (is not optimized); otherwise leave it unchanged

  @return returns 0 upon success and non-zero otherwise

  @see ComputeDotProduct_ref
*/
int ComputeDotProduct(const local_int_t start, const local_int_t end, const Vector & x, const Vector & y,
    double & result, bool & isOptimized) {

  (void) isOptimized;

  assert(x.localLength>=end); // Test vector lengths
  assert(y.localLength>=end);

  double local_result = 0.0;
  double * xv = x.values;
  double * yv = y.values;
  if (yv==xv) {
    for (local_int_t i=start; i<end; i++) local_result += xv[i]*xv[i];
  } else {
    for (local_int_t i=start; i<end; i++) local_result += xv[i]*yv[i];
  }

  result += local_result;

  return 0;
}
