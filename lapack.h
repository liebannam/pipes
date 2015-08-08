/* BLAS ROUTINES */

#ifndef JON_COMPLEX
#define JON_COMPLEX
#ifndef HIDE_FROM_DOXYGEN
typedef struct single_cx_struct { float r, i; } single_cx;
typedef struct double_cx_struct { double r, i; } double_cx;
#endif
#endif

#ifdef __cplusplus
#include <complex>
typedef std::complex<double> Complex;
extern "C" {
#else
typedef double_cx Complex;
#endif

/* dscal solves x = ax */
void dscal(int n, double a, double *x, int incx);
/* daxpy solves y = ax + y */
void  daxpy(int n, double a, double *x, int incx, double *y, int incy);
void  zaxpy(int n, Complex a, Complex *x, int incx, Complex *y, int incy);
/* dcopy replaces y with x */
void  dcopy(int n, double *x, int incx, double *y, int incy);
/* ddot returns the dot product of x and y */
double ddot(int n, double *x, int incx, double *y, int incy);
double dnrm2(int n, double *x, int incx);
int idamax(int n, double *x, int incx);
void dtpsv(char uplo, char transA, char diag, int n, double *A,
	   double *y, int incy);
/* matrix vector multiply */
void dgemv(char *trans, int *m, int *n, double *alpha, double *A, int *lda,
	   double *x, int *incx, double *beta, double *y, int *incy);
/* C = alpha op(A)*op(B) + beta C, op(A) MxK, op(B) KxN, op(X) = X or X' */
void dgemm(char transA, char transB, int m, int n, int k, double alpha,
	   double *A, int lda, double *B, int ldb, double beta, double *C,
	   int ldc);
/* C = alpha op(A)*op(B) + beta C, op(A) MxK, op(B) KxN, op(X) = X or X' */
void zgemm(char transA, char transB, int m, int n, int k, Complex alpha,
	   Complex *A, int lda, Complex *B, int ldb, Complex beta,
	   Complex *C, int ldc);
/* C = alpha AA' + beta C, alpha A'A + beta C, C is n x n */
void dsyrk(char uplo, char trans, int n, int k, double alpha,
	   double *A, int lda, double beta, double *C, int ldc);
/* B = X, the solution of op(A)*X = alpha*B  or  X*op(A) = alpha*B */
void dtrsm(char side, char uplo, char transA, char diag, int m, int n,
           double alpha, double *A, int ldA, double *B, int ldB);

/* LAPACK ROUTINES */
/* dgetrf replaces A by its LU decomposition */
void dgetrf(int m, int n, double *A, int lda, int *ipiv, int *info);
void zgetrf(int m, int n, Complex *A, int lda, int *ipiv, int *info);
/* dgetrs solves Ax=b; A should be LU decomposed, x replaces b on exit */
void dgetrs(char trans, int n, int nrhs, double *A, int lda,
	    int *ipiv, double *B, int ldb, int *info);
void zgetrs(char trans, int n, int nrhs, Complex *A, int lda,
	    int *ipiv, Complex *B, int ldb, int *info);
/* dgetri replaces the LU decomposition of A by its inverse */
void dgetri(int n, double *A, int lda, int *ipiv, int *info);
/* dgbtrf replaces A by its LU decomposition; A is banded */
void dgbtrf(int m, int n, int kl, int ku, double *ab, int ldab,
	    int *ipiv, int *info);
/* dgbtrs solves Ax=b with A banded and LU-decomposed */
void dgbtrs(char trans, int n, int kl, int ku, int nrhs, double *ab,
	    int ldab, int *ipiv, double *b, int ldb, int *info);
/* DGESVX uses the LU factorization to compute the solution to a real
 * system of linear equations A * X = B, where A is an N-by-N matrix
 * and X and B are N-by-NRHS matrices.  Error bounds on the solution
 * and a condition estimate are also provided. */
void dgesvx(char fact, char trans, int n, int nrhs, double *A, 
	    int lda, double *AF, int ldaf, int *ipiv, char *equed,
	    double *R, double *C, double *B, int ldb,
	    double *X, int ldx, double *rcond, double *ferr,
	    double *berr, int *info);

/* dgelsx computes the minimum norm least squares solution to an over- or
   under-determined system of linear equations AX=B, using a complete
   orthogonal factorization of A */
void dgelsx(int m, int n, int nrhs, double *A, int lda, double *B,
	    int ldb, int *jpvt, double rcond, int *rank, int *info);
void dgeqp3(int m, int n, double *A, int lda, int *jpvt,
	    double *tau, int *info);
void dgeqrf(int m, int n, double *A, int lda, double *tau, int *info);
void dormqr(char side, char trans, int M, int N, int K, double *Q, int ldq,
	    double *tau, double *C, int ldc, int *info);

/* dgelss compute a minimum-norm solution to a linear least squares
   problem minimize || B-AX ||2 using the SVD of a general matrix A. */
void dgelss(int m, int n, int nrhs, double *A, int lda, double *B,
	    int ldb, double *S, double rcond, int *rank, int *info);
void zgelss(int m, int n, int nrhs, Complex *A, int lda, Complex *B,
	    int ldb, double *S, double rcond, int *rank, int *info);
void dposv(char uplo, int n, int nrhs, double *A, int lda,
	   double *B, int ldb, int *info);
void dposvx(char fact, char uplo, int n, int nrhs, double *A, 
	    int lda, double *AF, int ldaf, char *equed,
	    double *S, double *B, int ldb,
	    double *X, int ldx, double *rcond, double *ferr,
	    double *berr, int *info);

void dpotrf(char uplo, int n, double *A, int lda, int *info);
void dpotrs(char uplo, int n, int nrhs,
	    double *A, int lda, double *B, int ldb, int *info);

void dgees(char jobv, char sort, int (*select)(double *),
	   int n, double *A, int lda, int *sdim, double *WR,
	   double *WI, double *VS, int ldvs, int *info);
void dgeev(char jobvl, char jobvr, int n, double *A, int lda,
	   double *WR, double *WI, double *VL, int ldvl,
	   double *VR, int ldvr, int *info);
void zgeev(char jobvl, char jobvr, int n, Complex *A, int lda,
	   Complex *W, Complex *VL, int ldvl,
	   Complex *VR, int ldvr, int *info);

void dgesvd(char* jobu, char *jobvt, int *m, int *n, double *A, int *lda,
	    double *S, double *U, int *ldu, double *VT, int *ldvt, 
	    double *work, int *lwork, int *info);
void zgesvd(char* jobu, char *jobvt, int *m, int *n, Complex *A, int *lda,
	    double *S, Complex *U, int *ldu, Complex *VT, int *ldvt,
	    double *work, int *lwork, int *info);

void dsyev(char jobz, char uplo, int n, double *A, int lda,
	   double *W, int *info);
void zheev(char jobz, char uplo, int n, Complex *A, int lda,
	   double *W, int *info);


#ifdef __cplusplus
}
#endif
