#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "strassen.h"


unsigned final_depth;

/***********************************************************************
 * Naive sequential algorithm, for comparison purposes
 **********************************************************************/
void matrixmul(int n, REAL *A, int an, REAL *B, int bn, REAL *C, int cn)
{
	int i, j, k;
	REAL s;

	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			s = 0.0;
			for (k = 0; k < n; ++k)
				s += ELEM(A, an, i, k) * ELEM(B, bn, k, j);

			ELEM(C, cn, i, j) = s;
		}
	}
}

/*****************************************************************************
**
** FastNaiveMatrixMultiply
**
** For small to medium sized matrices A, B, and C of size
** MatrixSize * MatrixSize this function performs the operation
** C = A x B efficiently.
**
** Note MatrixSize must be divisible by 8.
**
** INPUT:
**    C = (*C WRITE) Address of top left element of matrix C.
**    A = (*A IS READ ONLY) Address of top left element of matrix A.
**    B = (*B IS READ ONLY) Address of top left element of matrix B.
**    MatrixSize = Size of matrices (for n*n matrix, MatrixSize = n)
**    RowWidthA = Number of elements in memory between A[x,y] and A[x,y+1]
**    RowWidthB = Number of elements in memory between B[x,y] and B[x,y+1]
**    RowWidthC = Number of elements in memory between C[x,y] and C[x,y+1]
**
** OUTPUT:
**    C = (*C WRITE) Matrix C contains A x B. (Initial value of *C undefined.)
**
*****************************************************************************/
void FastNaiveMatrixMultiply(REAL *C, REAL *A, REAL *B, unsigned MatrixSize,
     unsigned RowWidthC, unsigned RowWidthA, unsigned RowWidthB)
{
	/* Assumes size of real is 8 bytes */
	PTR RowWidthBInBytes = RowWidthB  << 3;
	PTR RowWidthAInBytes = RowWidthA << 3;
	PTR MatrixWidthInBytes = MatrixSize << 3;
	PTR RowIncrementC = ( RowWidthC - MatrixSize) << 3;
	unsigned Horizontal, Vertical;

	REAL *ARowStart = A;
	for (Vertical = 0; Vertical < MatrixSize; Vertical++) {
		for (Horizontal = 0; Horizontal < MatrixSize; Horizontal += 8) {
			REAL *BColumnStart = B + Horizontal;
			REAL FirstARowValue = *ARowStart++;

			REAL Sum0 = FirstARowValue * (*BColumnStart);
			REAL Sum1 = FirstARowValue * (*(BColumnStart+1));
			REAL Sum2 = FirstARowValue * (*(BColumnStart+2));
			REAL Sum3 = FirstARowValue * (*(BColumnStart+3));
			REAL Sum4 = FirstARowValue * (*(BColumnStart+4));
			REAL Sum5 = FirstARowValue * (*(BColumnStart+5));
			REAL Sum6 = FirstARowValue * (*(BColumnStart+6));
			REAL Sum7 = FirstARowValue * (*(BColumnStart+7));

			unsigned Products;
			for (Products = 1; Products < MatrixSize; Products++) {
				REAL ARowValue = *ARowStart++;
				BColumnStart = (REAL*) (((PTR) BColumnStart) + RowWidthBInBytes);
				Sum0 += ARowValue * (*BColumnStart);
				Sum1 += ARowValue * (*(BColumnStart+1));
				Sum2 += ARowValue * (*(BColumnStart+2));
				Sum3 += ARowValue * (*(BColumnStart+3));
				Sum4 += ARowValue * (*(BColumnStart+4));
				Sum5 += ARowValue * (*(BColumnStart+5));
				Sum6 += ARowValue * (*(BColumnStart+6));
				Sum7 += ARowValue * (*(BColumnStart+7));
			}
			ARowStart = (REAL*) ( ((PTR) ARowStart) - MatrixWidthInBytes);

			*(C) = Sum0;
			*(C+1) = Sum1;
			*(C+2) = Sum2;
			*(C+3) = Sum3;
			*(C+4) = Sum4;
			*(C+5) = Sum5;
			*(C+6) = Sum6;
			*(C+7) = Sum7;
			C+=8;
		}
		ARowStart = (REAL*) ( ((PTR) ARowStart) + RowWidthAInBytes );
		C = (REAL*) ( ((PTR) C) + RowIncrementC );
	}
}

/*****************************************************************************
**
** FastAdditiveNaiveMatrixMultiply
**
** For small to medium sized matrices A, B, and C of size
** MatrixSize * MatrixSize this function performs the operation
** C += A x B efficiently.
**
** Note MatrixSize must be divisible by 8.
**
** INPUT:
**    C = (*C READ/WRITE) Address of top left element of matrix C.
**    A = (*A IS READ ONLY) Address of top left element of matrix A.
**    B = (*B IS READ ONLY) Address of top left element of matrix B.
**    MatrixSize = Size of matrices (for n*n matrix, MatrixSize = n)
**    RowWidthA = Number of elements in memory between A[x,y] and A[x,y+1]
**    RowWidthB = Number of elements in memory between B[x,y] and B[x,y+1]
**    RowWidthC = Number of elements in memory between C[x,y] and C[x,y+1]
**
** OUTPUT:
**    C = (*C READ/WRITE) Matrix C contains C + A x B.
**
*****************************************************************************/
void FastAdditiveNaiveMatrixMultiply(REAL *C, REAL *A, REAL *B, unsigned MatrixSize,
     unsigned RowWidthC, unsigned RowWidthA, unsigned RowWidthB)
{
	/* Assumes size of real is 8 bytes */
	PTR RowWidthBInBytes = RowWidthB  << 3;
	PTR RowWidthAInBytes = RowWidthA << 3;
	PTR MatrixWidthInBytes = MatrixSize << 3;
	PTR RowIncrementC = ( RowWidthC - MatrixSize) << 3;
	unsigned Horizontal, Vertical;

	REAL *ARowStart = A;
	for (Vertical = 0; Vertical < MatrixSize; Vertical++) {
		for (Horizontal = 0; Horizontal < MatrixSize; Horizontal += 8) {
			REAL *BColumnStart = B + Horizontal;

			REAL Sum0 = *C;
			REAL Sum1 = *(C+1);
			REAL Sum2 = *(C+2);
			REAL Sum3 = *(C+3);
			REAL Sum4 = *(C+4);
			REAL Sum5 = *(C+5);
			REAL Sum6 = *(C+6);
			REAL Sum7 = *(C+7);

			unsigned Products;
			for (Products = 0; Products < MatrixSize; Products++) {
				REAL ARowValue = *ARowStart++;

				Sum0 += ARowValue * (*BColumnStart);
				Sum1 += ARowValue * (*(BColumnStart+1));
				Sum2 += ARowValue * (*(BColumnStart+2));
				Sum3 += ARowValue * (*(BColumnStart+3));
				Sum4 += ARowValue * (*(BColumnStart+4));
				Sum5 += ARowValue * (*(BColumnStart+5));
				Sum6 += ARowValue * (*(BColumnStart+6));
				Sum7 += ARowValue * (*(BColumnStart+7));

				BColumnStart = (REAL*) (((PTR) BColumnStart) + RowWidthBInBytes);
			}
			ARowStart = (REAL*) ( ((PTR) ARowStart) - MatrixWidthInBytes);

			*(C) = Sum0;
			*(C+1) = Sum1;
			*(C+2) = Sum2;
			*(C+3) = Sum3;
			*(C+4) = Sum4;
			*(C+5) = Sum5;
			*(C+6) = Sum6;
			*(C+7) = Sum7;
			C+=8;
		}

		ARowStart = (REAL*) ( ((PTR) ARowStart) + RowWidthAInBytes );
		C = (REAL*) ( ((PTR) C) + RowIncrementC );
	}
}

/*****************************************************************************
**
** MultiplyByDivideAndConquer
**
** For medium to medium-large (would you like fries with that) sized
** matrices A, B, and C of size MatrixSize * MatrixSize this function
** efficiently performs the operation
**    C  = A x B (if AdditiveMode == 0)
**    C += A x B (if AdditiveMode != 0)
**
** Note MatrixSize must be divisible by 16.
**
** INPUT:
**    C = (*C READ/WRITE) Address of top left element of matrix C.
**    A = (*A IS READ ONLY) Address of top left element of matrix A.
**    B = (*B IS READ ONLY) Address of top left element of matrix B.
**    MatrixSize = Size of matrices (for n*n matrix, MatrixSize = n)
**    RowWidthA = Number of elements in memory between A[x,y] and A[x,y+1]
**    RowWidthB = Number of elements in memory between B[x,y] and B[x,y+1]
**    RowWidthC = Number of elements in memory between C[x,y] and C[x,y+1]
**    AdditiveMode = 0 if we want C = A x B, otherwise we'll do C += A x B
**
** OUTPUT:
**    C (+)= A x B. (+ if AdditiveMode != 0)
**
*****************************************************************************/
void MultiplyByDivideAndConquer(REAL *C, REAL *A, REAL *B,
	unsigned MatrixSize,
	unsigned RowWidthC,
	unsigned RowWidthA,
	unsigned RowWidthB,
	int AdditiveMode
)
{
	#define A00 A
	#define B00 B
	#define C00 C
	REAL  *A01, *A10, *A11, *B01, *B10, *B11, *C01, *C10, *C11;
	unsigned QuadrantSize = MatrixSize >> 1;

	/* partition the matrix */
	A01 = A00 + QuadrantSize;
	A10 = A00 + RowWidthA * QuadrantSize;
	A11 = A10 + QuadrantSize;

	B01 = B00 + QuadrantSize;
	B10 = B00 + RowWidthB * QuadrantSize;
	B11 = B10 + QuadrantSize;

	C01 = C00 + QuadrantSize;
	C10 = C00 + RowWidthC * QuadrantSize;
	C11 = C10 + QuadrantSize;

	if (QuadrantSize > SizeAtWhichNaiveAlgorithmIsMoreEfficient)
	{
		MultiplyByDivideAndConquer(C00, A00, B00, QuadrantSize,
			RowWidthC, RowWidthA, RowWidthB, AdditiveMode);

		MultiplyByDivideAndConquer(C01, A00, B01, QuadrantSize,
			RowWidthC, RowWidthA, RowWidthB, AdditiveMode);

		MultiplyByDivideAndConquer(C11, A10, B01, QuadrantSize,
			RowWidthC, RowWidthA, RowWidthB, AdditiveMode);

		MultiplyByDivideAndConquer(C10, A10, B00, QuadrantSize,
			RowWidthC, RowWidthA, RowWidthB, AdditiveMode);

		MultiplyByDivideAndConquer(C00, A01, B10, QuadrantSize,
			RowWidthC, RowWidthA, RowWidthB, 1);

		MultiplyByDivideAndConquer(C01, A01, B11, QuadrantSize,
			RowWidthC, RowWidthA, RowWidthB, 1);

		MultiplyByDivideAndConquer(C11, A11, B11, QuadrantSize,
			RowWidthC, RowWidthA, RowWidthB, 1);

		MultiplyByDivideAndConquer(C10, A11, B10, QuadrantSize,
			RowWidthC, RowWidthA, RowWidthB, 1);
	}
	else {
		if (AdditiveMode) {
			FastAdditiveNaiveMatrixMultiply(C00, A00, B00, QuadrantSize,
				RowWidthC, RowWidthA, RowWidthB);

			FastAdditiveNaiveMatrixMultiply(C01, A00, B01, QuadrantSize,
				RowWidthC, RowWidthA, RowWidthB);

			FastAdditiveNaiveMatrixMultiply(C11, A10, B01, QuadrantSize,
				RowWidthC, RowWidthA, RowWidthB);

			FastAdditiveNaiveMatrixMultiply(C10, A10, B00, QuadrantSize,
				RowWidthC, RowWidthA, RowWidthB);
		}
		else {
			FastNaiveMatrixMultiply(C00, A00, B00, QuadrantSize,
				RowWidthC, RowWidthA, RowWidthB);

			FastNaiveMatrixMultiply(C01, A00, B01, QuadrantSize,
				RowWidthC, RowWidthA, RowWidthB);

			FastNaiveMatrixMultiply(C11, A10, B01, QuadrantSize,
				RowWidthC, RowWidthA, RowWidthB);

			FastNaiveMatrixMultiply(C10, A10, B00, QuadrantSize,
				RowWidthC, RowWidthA, RowWidthB);
		}

		FastAdditiveNaiveMatrixMultiply(C00, A01, B10, QuadrantSize,
			RowWidthC, RowWidthA, RowWidthB);

		FastAdditiveNaiveMatrixMultiply(C01, A01, B11, QuadrantSize,
			RowWidthC, RowWidthA, RowWidthB);

		FastAdditiveNaiveMatrixMultiply(C11, A11, B11, QuadrantSize,
			RowWidthC, RowWidthA, RowWidthB);

		FastAdditiveNaiveMatrixMultiply(C10, A11, B10, QuadrantSize,
			RowWidthC, RowWidthA, RowWidthB);
	}
	return;
}

void OptimizedStrassenMultiply_par(REAL *C, REAL *A, REAL *B, unsigned MatrixSize,
     unsigned RowWidthC, unsigned RowWidthA, unsigned RowWidthB, int Depth)
{
	unsigned QuadrantSize = MatrixSize >> 1; /* MatrixSize / 2 */
	unsigned QuadrantSizeInBytes = sizeof(REAL) * QuadrantSize * QuadrantSize + 32;
	unsigned Column, Row;

	/************************************************************************
	** For each matrix A, B, and C, we'll want pointers to each quandrant
	** in the matrix. These quandrants will be addressed as follows:
	**  --        --
	**  | A11  A12 |
	**  |          |
	**  | A21  A22 |
	**  --        --
	************************************************************************/
	REAL /* *A11, *B11, *C11, */ *A12, *B12, *C12, *A21, *B21, *C21, *A22, *B22, *C22;
	REAL *S1,*S2,*S3,*S4,*S5,*S6,*S7,*S8,*M2,*M5,*T1sMULT;
	#define T2sMULT C22
	#define NumberOfVariables 11

	PTR TempMatrixOffset = 0;
	PTR MatrixOffsetA = 0;
	PTR MatrixOffsetB = 0;

	char *Heap;
	void *StartHeap;

	/* Distance between the end of a matrix row and the start of the next row */
	PTR RowIncrementA = ( RowWidthA - QuadrantSize ) << 3;
	PTR RowIncrementB = ( RowWidthB - QuadrantSize ) << 3;
	PTR RowIncrementC = ( RowWidthC - QuadrantSize ) << 3;

	if (MatrixSize <= 64) {
		MultiplyByDivideAndConquer(C, A, B, MatrixSize, RowWidthC, RowWidthA, RowWidthB, 0);
		return;
	}

	/* Initialize quadrant matrices */
	#define A11 A
	#define B11 B
	#define C11 C
	A12 = A11 + QuadrantSize;
	B12 = B11 + QuadrantSize;
	C12 = C11 + QuadrantSize;
	A21 = A + (RowWidthA * QuadrantSize);
	B21 = B + (RowWidthB * QuadrantSize);
	C21 = C + (RowWidthC * QuadrantSize);
	A22 = A21 + QuadrantSize;
	B22 = B21 + QuadrantSize;
	C22 = C21 + QuadrantSize;

	/* Allocate Heap Space Here */
	StartHeap = Heap = malloc(QuadrantSizeInBytes * NumberOfVariables);
	/* ensure that heap is on cache boundary */
	if ( ((PTR) Heap) & 31)
		Heap = (char*) ( ((PTR) Heap) + 32 - ( ((PTR) Heap) & 31) );

	/* Distribute the heap space over the variables */
	S1 = (REAL*) Heap; Heap += QuadrantSizeInBytes;
	S2 = (REAL*) Heap; Heap += QuadrantSizeInBytes;
	S3 = (REAL*) Heap; Heap += QuadrantSizeInBytes;
	S4 = (REAL*) Heap; Heap += QuadrantSizeInBytes;
	S5 = (REAL*) Heap; Heap += QuadrantSizeInBytes;
	S6 = (REAL*) Heap; Heap += QuadrantSizeInBytes;
	S7 = (REAL*) Heap; Heap += QuadrantSizeInBytes;
	S8 = (REAL*) Heap; Heap += QuadrantSizeInBytes;
	M2 = (REAL*) Heap; Heap += QuadrantSizeInBytes;
	M5 = (REAL*) Heap; Heap += QuadrantSizeInBytes;
	T1sMULT = (REAL*) Heap; Heap += QuadrantSizeInBytes;

	/***************************************************************************
	** Step through all columns row by row (vertically)
	** (jumps in memory by RowWidth => bad locality)
	** (but we want the best locality on the innermost loop)
	***************************************************************************/
	for (Row = 0; Row < QuadrantSize; Row++) {

		/*************************************************************************
		** Step through each row horizontally (addressing elements in each column)
		** (jumps linearly througn memory => good locality)
		*************************************************************************/
		for (Column = 0; Column < QuadrantSize; Column++) {

			/***********************************************************
			** Within this loop, the following holds for MatrixOffset:
			** MatrixOffset = (Row * RowWidth) + Column
			** (note: that the unit of the offset is number of reals)
			***********************************************************/
			/* Element of Global Matrix, such as A, B, C */
			#define E(Matrix)   (* (REAL*) ( ((PTR) Matrix) + TempMatrixOffset ) )
			#define EA(Matrix)  (* (REAL*) ( ((PTR) Matrix) + MatrixOffsetA ) )
			#define EB(Matrix)  (* (REAL*) ( ((PTR) Matrix) + MatrixOffsetB ) )

			/* FIXME - may pay to expand these out - got higher speed-ups below */
			/* S4 = A12 - ( S2 = ( S1 = A21 + A22 ) - A11 ) */
			E(S4) = EA(A12) - ( E(S2) = ( E(S1) = EA(A21) + EA(A22) ) - EA(A11) );

			/* S8 = (S6 = B22 - ( S5 = B12 - B11 ) ) - B21 */
			E(S8) = ( E(S6) = EB(B22) - ( E(S5) = EB(B12) - EB(B11) ) ) - EB(B21);

			/* S3 = A11 - A21 */
			E(S3) = EA(A11) - EA(A21);

			/* S7 = B22 - B12 */
			E(S7) = EB(B22) - EB(B12);

			TempMatrixOffset += sizeof(REAL);
			MatrixOffsetA += sizeof(REAL);
			MatrixOffsetB += sizeof(REAL);
		} /* end row loop*/

		MatrixOffsetA += RowIncrementA;
		MatrixOffsetB += RowIncrementB;
	} /* end column loop */

	unsigned long long algo_cost = (unsigned long long) pow(QuadrantSize, 2.8074);

	/* M2 = A11 x B11 */
	#pragma oss task label("strassen") cost(algo_cost) final(QuadrantSize <= final_depth)
	OptimizedStrassenMultiply_par(M2, A11, B11, QuadrantSize, QuadrantSize, RowWidthA, RowWidthB, Depth+1);

	/* M5 = S1 * S5 */
	#pragma oss task label("strassen") cost(algo_cost) final(QuadrantSize <= final_depth)
	OptimizedStrassenMultiply_par(M5, S1, S5, QuadrantSize, QuadrantSize, QuadrantSize, QuadrantSize, Depth+1);

	/* Step 1 of T1 = S2 x S6 + M2 */
	#pragma oss task label("strassen") cost(algo_cost) final(QuadrantSize <= final_depth)
	OptimizedStrassenMultiply_par(T1sMULT, S2, S6,  QuadrantSize, QuadrantSize, QuadrantSize, QuadrantSize, Depth+1);

	/* Step 1 of T2 = T1 + S3 x S7 */
	#pragma oss task label("strassen") cost(algo_cost) final(QuadrantSize <= final_depth)
	OptimizedStrassenMultiply_par(C22, S3, S7, QuadrantSize, RowWidthC /*FIXME*/, QuadrantSize, QuadrantSize, Depth+1);

	/* Step 1 of C11 = M2 + A12 * B21 */
	#pragma oss task label("strassen") cost(algo_cost) final(QuadrantSize <= final_depth)
	OptimizedStrassenMultiply_par(C11, A12, B21, QuadrantSize, RowWidthC, RowWidthA, RowWidthB, Depth+1);

	/* Step 1 of C12 = S4 x B22 + T1 + M5 */
	#pragma oss task label("strassen") cost(algo_cost) final(QuadrantSize <= final_depth)
	OptimizedStrassenMultiply_par(C12, S4, B22, QuadrantSize, RowWidthC, QuadrantSize, RowWidthB, Depth+1);

	/* Step 1 of C21 = T2 - A22 * S8 */
	#pragma oss task label("strassen") cost(algo_cost) final(QuadrantSize <= final_depth)
	OptimizedStrassenMultiply_par(C21, A22, S8, QuadrantSize, RowWidthC, RowWidthA, QuadrantSize, Depth+1);

	/**********************************************
	** Synchronization Point
	**********************************************/
	#pragma oss taskwait
	/***************************************************************************
	** Step through all columns row by row (vertically)
	** (jumps in memory by RowWidth => bad locality)
	** (but we want the best locality on the innermost loop)
	***************************************************************************/
	for (Row = 0; Row < QuadrantSize; Row++) {
		/*************************************************************************
		** Step through each row horizontally (addressing elements in each column)
		** (jumps linearly througn memory => good locality)
		*************************************************************************/
		for (Column = 0; Column < QuadrantSize; Column += 4) {
			REAL LocalM5_0 = *(M5  );
			REAL LocalM5_1 = *(M5+1);
			REAL LocalM5_2 = *(M5+2);
			REAL LocalM5_3 = *(M5+3);
			REAL LocalM2_0 = *(M2  );
			REAL LocalM2_1 = *(M2+1);
			REAL LocalM2_2 = *(M2+2);
			REAL LocalM2_3 = *(M2+3);
			REAL T1_0 = *(T1sMULT  ) + LocalM2_0;
			REAL T1_1 = *(T1sMULT+1) + LocalM2_1;
			REAL T1_2 = *(T1sMULT+2) + LocalM2_2;
			REAL T1_3 = *(T1sMULT+3) + LocalM2_3;
			REAL T2_0 = *(C22  ) + T1_0;
			REAL T2_1 = *(C22+1) + T1_1;
			REAL T2_2 = *(C22+2) + T1_2;
			REAL T2_3 = *(C22+3) + T1_3;
			(*(C11  )) += LocalM2_0;
			(*(C11+1)) += LocalM2_1;
			(*(C11+2)) += LocalM2_2;
			(*(C11+3)) += LocalM2_3;
			(*(C12  )) += LocalM5_0 + T1_0;
			(*(C12+1)) += LocalM5_1 + T1_1;
			(*(C12+2)) += LocalM5_2 + T1_2;
			(*(C12+3)) += LocalM5_3 + T1_3;
			(*(C22  ))  = LocalM5_0 + T2_0;
			(*(C22+1))  = LocalM5_1 + T2_1;
			(*(C22+2))  = LocalM5_2 + T2_2;
			(*(C22+3))  = LocalM5_3 + T2_3;
			(*(C21  ))  = (- *(C21  )) + T2_0;
			(*(C21+1))  = (- *(C21+1)) + T2_1;
			(*(C21+2))  = (- *(C21+2)) + T2_2;
			(*(C21+3))  = (- *(C21+3)) + T2_3;
			M5 += 4;
			M2 += 4;
			T1sMULT += 4;
			C11 += 4;
			C12 += 4;
			C21 += 4;
			C22 += 4;
		}
		C11 = (REAL*) ( ((PTR) C11 ) + RowIncrementC);
		C12 = (REAL*) ( ((PTR) C12 ) + RowIncrementC);
		C21 = (REAL*) ( ((PTR) C21 ) + RowIncrementC);
		C22 = (REAL*) ( ((PTR) C22 ) + RowIncrementC);
	}

	free(StartHeap);
}

void init_matrix(int n, REAL *A, int an)
{
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			ELEM(A, an, i, j) = ((double) rand()) / (double) RAND_MAX;
		}
	}
}


int main(int argc, char **argv)
{
	// Initialize
	int reps    = 3;
	unsigned n  = 2048;
	final_depth = 1024;

	if (argc > 4 || (argc >= 2 && strcmp(argv[1], "-h") == 0)) {
		printf("\nUsage: ./strassen [-h] <R> <N> <F>\n");
		printf("  R : Number of repetitions\n");
		printf("  N : Matrix dimension (N*N). MUST be a power of 2 and multiple of 16\n");
		printf("  F : The maximum depth in which tasks are final (chunksize)\n\n");
		exit(1);
	}

	if (argc >= 2)
		reps = atoi(argv[1]);
	if (argc >= 3)
		n = atoi(argv[2]);
	if (argc >= 4)
		final_depth = atoi(argv[3]);

	if ((n & (n - 1)) != 0 || (n % 16) != 0) {
		printf("Error: matrix size (%d) must be a power of 2 and a multiple of %d\n", n, 16);
		exit (1);
	}

	double * A = (double *) malloc(n * n * sizeof(double));
	double * B = (double *) malloc(n * n * sizeof(double));
	double * C = (double *) malloc(n * n * sizeof(double));

	init_matrix(n, A, n);
	init_matrix(n, B, n);

	// Warmup iteration
	#pragma oss task label("strassen") cost((unsigned long long) pow(n, 2.8074)) final(n <= final_depth)
	OptimizedStrassenMultiply_par(C, A, B, n, n, n, n, 1);
	#pragma oss taskwait

	// Start timing
	struct timeval start;
	gettimeofday(&start, NULL);

	// Main loop
	for (int i = 0; i < reps; ++i) {
		#pragma oss task label("strassen") cost((unsigned long long) pow(n, 2.8074)) final(n <= final_depth)
		OptimizedStrassenMultiply_par(C, A, B, n, n, n, n, 1);
		#pragma oss taskwait
	}

	// Stop timing
	struct timeval stop;
	gettimeofday(&stop, NULL);
	unsigned long elapsed = 1000000 * (stop.tv_sec - start.tv_sec);
	elapsed += stop.tv_usec - start.tv_usec;
	double t = (double) elapsed * 1e-6;
	double perf = ((pow(n, 2.8074)) / ((double)elapsed/1000000.0))/1000000.0;

	// Print statistics
	printf("%14e %14e %14d %14d %14d\n", t, perf, reps, n, final_depth);

	// Destroy
	free(A);
	free(B);
	free(C);
}

