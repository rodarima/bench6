
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

#ifndef COMPUTESYMGS_HPP
#define COMPUTESYMGS_HPP
#include "SparseMatrix.hpp"
#include "Vector.hpp"

int ComputeSYMGS( const SparseMatrix & A, const Vector & r, Vector & x);
int ComputeFGS( const SparseMatrix & A, local_int_t start, local_int_t end, const Vector & r, Vector & x);
int ComputeBGS( const SparseMatrix & A, local_int_t start, local_int_t end, const Vector & r, Vector & x);

#endif // COMPUTESYMGS_HPP
