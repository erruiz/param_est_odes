#ifndef BLAS_H
#define BLAS_H

#define F77_CALL(x) x ## _
#define F77_NAME(x) F77_CALL(x)

#ifdef __cplusplus
extern "C" {
#endif

/* Level 1 BLAS */
extern void F77_NAME(dscal)(const int *n,const double *da,const double *dx,
	const int *incx);

/* Level 2 BLAS */
extern void F77_NAME(dgemv)(const char *trans,const int *m,const int *n,
	const double *alpha,const double *a,const int *lda,const double *x,
	const int *incx,const double *beta,const double *y,const int *incy);

/* Level 3 BLAS */
extern void F77_NAME(dgemm)(const char *transa,const char *transb,
	const int *m, const int *n, const int *k,
	const double *alpha, const double *a, const int *lda,const double *b, 
	const int *ldb, const double *beta, double *c, const int *ldc);

#ifdef __cplusplus
}
#endif

#endif
