
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

#ifndef EXCHANGEHALOREF_HPP
#define EXCHANGEHALOREF_HPP
#include "SparseMatrix.hpp"
#include "Vector.hpp"
void ExchangeHalo_ref(const SparseMatrix & A, Vector & x);
#endif // EXCHANGEHALO_HPP
