
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

#ifndef COMPUTERESTRICTION_HPP
#define COMPUTERESTRICTION_HPP
#include "Vector.hpp"
#include "SparseMatrix.hpp"
int ComputeRestriction(const SparseMatrix & A, local_int_t start, local_int_t end, const Vector & rf);
#endif // COMPUTERESTRICTION_REF_HPP
