#ifndef _STRASSEN_H
#define _STRASSEN_H

typedef double REAL;
typedef unsigned long PTR;

/* FIXME: At the moment we use a constant value, change to parameter ???*/
/* Below this cut off strassen uses FastAdditiveNaiveMatrixMultiply algorithm */
#define SizeAtWhichNaiveAlgorithmIsMoreEfficient 16

/***********************************************************************
 * Maximum tolerable relative error (for the checking routine)
 **********************************************************************/
#define EPSILON (1.0E-6)
/***********************************************************************
 * Matrices are stored in row-major order; A is a pointer to
 * the first element of the matrix, and an is the number of elements
 * between two rows. This macro produces the element A[i,j]
 * given A, an, i and j
 **********************************************************************/
#define ELEM(A, an, i, j) (A[(i)*(an)+(j)])

void matrixmul(int n, REAL *A, int an, REAL *B, int bn, REAL *C, int cn);

void FastNaiveMatrixMultiply(REAL *C, REAL *A, REAL *B, unsigned MatrixSize,
	unsigned RowWidthC, unsigned RowWidthA, unsigned RowWidthB);

void FastAdditiveNaiveMatrixMultiply(REAL *C, REAL *A, REAL *B, unsigned MatrixSize,
	unsigned RowWidthC, unsigned RowWidthA, unsigned RowWidthB);

void MultiplyByDivideAndConquer(REAL *C, REAL *A, REAL *B,
	unsigned MatrixSize,
	unsigned RowWidthC,
	unsigned RowWidthA,
	unsigned RowWidthB,
	int AdditiveMode
);

void OptimizedStrassenMultiply_par(REAL *C, REAL *A, REAL *B, unsigned MatrixSize,
	unsigned RowWidthC, unsigned RowWidthA, unsigned RowWidthB, int Depth);

void init_matrix(int n, REAL *A, int an);

#endif // _STRASSEN_H

