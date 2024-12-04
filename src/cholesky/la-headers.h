#if HAVE_CBLAS_H
#include <cblas.h>
#elif HAVE_MKL_CBLAS_H
#include <mkl_cblas.h>
#else
#error Need a CBLAS library
#endif

#if HAVE_LAPACKE_H
#include <lapacke.h>

#elif HAVE_MKL_LAPACKE_H
#include <mkl_lapacke.h>

#elif HAVE_CLAPACK_H
#include <clapack.h>
#include <stdlib.h>
	/* Hand-made wrapper */
	typedef int lapack_int;
	#define LAPACK_ROW_MAJOR ((int) CblasRowMajor)
	
	static inline lapack_int LAPACKE_dpotrf(int matrix_layout, char uplo, lapack_int n, double *a, lapack_int lda)
	{
		enum ATLAS_UPLO clapack_uplo;

		if (uplo == 'L' || uplo == 'l') {
			clapack_uplo = CblasLower;
		} else if (uplo == 'U' || uplo == 'u') {
			clapack_uplo = CblasUpper;
		} else {
			abort();
		}
		return clapack_dpotrf((enum ATLAS_ORDER) matrix_layout, clapack_uplo, n, a, lda);
	}

#else
#error Need either LAPACKE or CLAPACK
#endif

