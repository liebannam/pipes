// LAPACK STUFF ////////////////

#include <stdio.h>
#include <stdlib.h>
#include "lapack.h"
#include "f_blas_lapack.h"

#define MAX(X,Y)  ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y)  ((X) < (Y) ? (X) : (Y))

// BLAS /////////////////////////

void dscal(int n, double a, double *x, int incx) {
  dscal_(&n, &a, x, &incx);
}
void  daxpy(int n, double a, double *x, int incx, double *y, int incy) {
  daxpy_(&n, &a, x, &incx, y, &incy);
}
void  zaxpy(int n, Complex a, Complex *x, int incx, Complex *y, int incy) {
  zaxpy_(&n, &a, x, &incx, y, &incy);
}
void  dcopy(int n, double *x, int incx, double *y, int incy) {
  dcopy_(&n, x, &incx, y, &incy);
}
double ddot(int n, double *x, int incx, double *y, int incy) {
  return ddot_(&n, x, &incx, y, &incy);
}
double dnrm2(int n, double *x, int incx) {
  return dnrm2_(&n, x, &incx);
}
int idamax(int n, double *x, int incx) {
  return idamax_(&n, x, &incx);
}
void dtpsv(char uplo, char transA, char diag, int n, double *A,
	   double *y, int incy) {
    dtpsv_(&uplo, &transA, &diag, &n, A, y, &incy);
}
void dgemv(char trans, int m, int n, double alpha, double *A, int lda,
	   double *x, int incx, double beta, double *y, int incy) {
  dgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}
void dgemm(char transA, char transB, int m, int n, int k, double alpha,
	   double *A, int lda, double *B, int ldb, double beta, double *C,
	   int ldc) {
  dgemm_(&transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb,
	 &beta, C, &ldc);
}
void zgemm(char transA, char transB, int m, int n, int k, Complex alpha,
	   Complex *A, int lda, Complex *B, int ldb, Complex beta,
	   Complex *C, int ldc) {
  zgemm_(&transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb,
	 &beta, C, &ldc);
}
void dsyrk(char uplo, char trans, int n, int k, double alpha,
	   double *A, int lda, double beta, double *C, int ldc) {
  dsyrk_(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc);
}

// LAPACK /////////////////////////////

void dgetrf(int m, int n, double *A, int lda, int *ipiv, int *info) {
  dgetrf_(&m, &n, A, &lda, ipiv, info);
}
void zgetrf(int m, int n, Complex *A, int lda, int *ipiv, int *info) {
  zgetrf_(&m, &n, A, &lda, ipiv, info);
}
void dgetrs(char trans, int n, int nrhs, double *A, int lda,
	    int *ipiv, double *B, int ldb, int *info) {
  dgetrs_(&trans, &n, &nrhs, A, &lda, ipiv, B, &ldb, info);
}
void zgetrs(char trans, int n, int nrhs, Complex *A, int lda,
	    int *ipiv, Complex *B, int ldb, int *info) {
  zgetrs_(&trans, &n, &nrhs, A, &lda, ipiv, B, &ldb, info);
}
void dgetri(int n, double *A, int lda, int *ipiv, int *info) {
  /* dim(WORK) = LWORK >= N */
  double *work = (double *) malloc (n*sizeof(double));
  dgetri_(&n, A, &lda, ipiv, work, &n, info);
  free(work);
}
void dgbtrf(int m, int n, int kl, int ku, double *ab, int ldab,
	    int *ipiv, int *info) {
  dgbtrf_(&m, &n, &kl, &ku, ab, &ldab, ipiv, info);
}
void dgbtrs(char trans, int n, int kl, int ku, int nrhs, double *ab,
	    int ldab, int *ipiv, double *b, int ldb, int *info) {
  dgbtrs_(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, info);
}
int dtrsm(char* side, char* uplo, char* trans, char* diag, int* M, int* N,
	  double* alpha, double* A, int* lda, double* B, int* ldb) {
	return dtrsm_(side,uplo,trans,diag,M,N,alpha,A,lda,B,ldb);
}
void dgesvx(char fact, char trans, int n, int nrhs, double *A, 
	    int lda, double *AF, int ldaf, int *ipiv, char *equed,
	    double *R, double *C, double *B, int ldb,
	    double *X, int ldx, double *rcond, double *ferr,
	    double *berr, int *info) {
  double *work = (double *) malloc (4*n * sizeof(double));
  int *iwork = (int *) malloc (n * sizeof(int));
  dgesvx_(&fact, &trans, &n, &nrhs, A, &lda, AF, &ldaf, ipiv,
	  equed, R, C, B, &ldb, X, &ldx, rcond, ferr, berr,
	  work, iwork, info);
  free(iwork);
  free(work);
}
void dgelsx(int m, int n, int nrhs, double *A, int lda, double *B,
	   int ldb, int *jpvt, double rcond, int *rank, int *info) {
  int worksize = MAX( MIN(m,n) + 3*n, 2*MIN(m,n)+nrhs ); 
  double *work = (double *) malloc (worksize*sizeof(double));
  dgelsx_(&m, &n, &nrhs, A, &lda, B, &ldb, jpvt, &rcond, rank, work, info);
  free(work);
}
void dgeqrf(int m, int n, double *A, int lda, double *tau, int *info)
{
    int lwork = -1;
    double *work = (double *) malloc (sizeof(double));
    dgeqrf_(&m, &n, A, &lda, tau, work, &lwork, info);
    lwork = (int) work[0];
    // printf("lwork = %d\n", lwork);
    work = (double *) realloc (work, lwork * sizeof(double));
    dgeqrf_(&m, &n, A, &lda, tau, work, &lwork, info);
    free(work);
}
void dormqr(char side, char trans, int M, int N, int K, double *Q, int ldq,
	    double *tau, double *C, int ldc, int *info) {
    int lwork = -1;
    double *work = (double *) malloc (sizeof(double));
    dormqr_(&side, &trans, &M, &N, &K, Q, &ldq, tau, C,
	    &ldc, work, &lwork, info);
    lwork = (int) work[0];
    // printf("lwork = %d\n", lwork);
    work = (double *) realloc (work, lwork * sizeof(double));
    dormqr_(&side, &trans, &M, &N, &K, Q, &ldq, tau, C,
	    &ldc, work, &lwork, info);
    free(work);
}
void dgelss(int m, int n, int nrhs, double *A, int lda, double *B,
	    int ldb, double *S, double rcond, int *rank,  int *info) {
    int lwork = -1;
    double *work = (double *) malloc (sizeof(double));
    if (nrhs<=0) {
	fprintf(stderr, "call zgelss with nrhs>0 (LAPACK bug)");
	exit(1);
    }
    dgelss_(&m, &n, &nrhs, A, &lda, B, &ldb, S, &rcond, rank,
	    work, &lwork, info);
    lwork = (int) work[0];
    work = (double *) realloc (work, lwork * sizeof(double));
    dgelss_(&m, &n, &nrhs, A, &lda, B, &ldb, S, &rcond, rank,
	    work, &lwork, info);
    free(work);
}
void zgelss(int m, int n, int nrhs, Complex *A, int lda, Complex *B,
	    int ldb, double *S, double rcond, int *rank,  int *info) {
    int lwork = -1;
    int rworksize = (m<n) ? 5*m : 5*n;
    Complex *work = (Complex *) malloc (sizeof(Complex));
    double *rwork = (double *) malloc (rworksize * sizeof(double));
    if (nrhs<=0) {
	fprintf(stderr, "call zgelss with nrhs>0 (LAPACK bug)\n");
	exit(1);
    }
    zgelss_(&m, &n, &nrhs, A, &lda, B, &ldb, S, &rcond, rank,
	    work, &lwork, rwork, info);
    lwork = (int) *((double *) work);
    //fprintf(stderr,"m: %d, n: %d, nrhs: %d, lwork: %d\n", m,n,nrhs,lwork);
    work = (Complex *) realloc (work, lwork * sizeof(Complex));
    zgelss_(&m, &n, &nrhs, A, &lda, B, &ldb, S, &rcond, rank,
	    work, &lwork, rwork, info);
    free(work);
    free(rwork);
}
void dposvx(char fact, char uplo, int n, int nrhs, double *A, 
	    int lda, double *AF, int ldaf, char *equed,
	    double *S, double *B, int ldb,
	    double *X, int ldx, double *rcond, double *ferr,
	    double *berr, int *info) {
  double *work = (double *) malloc (3*n * sizeof(double));
  int *iwork = (int *) malloc (n * sizeof(int));
  dposvx_(&fact, &uplo, &n, &nrhs, A, &lda, AF, &ldaf,
	  equed, S, B, &ldb, X, &ldx, rcond, ferr, berr,
	  work, iwork, info);
  free(iwork);
  free(work);
}
void dpotrf(char uplo, int n, double *A, int lda, int *info) {
    dpotrf_(&uplo, &n, A, &lda, info);
}
void dpotrs(char uplo, int n, int nrhs,
	    double *A, int lda, double *B, int ldb, int *info) {
    dpotrs_(&uplo, &n, &nrhs, A, &lda, B, &ldb, info);
}
void dgees2(char jobv, char sort, int (*select)(double *),
	   int n, double *A, int lda, int *sdim, double *WR,
	   double *WI, double *VS, int ldvs, int *info) {
    double *work = (double *) malloc (sizeof(double));
    int *bwork = (int *) malloc (n * sizeof(int));
    int lwork = -1;
    dgees_(&jobv, &sort, select, &n, A, &lda, sdim, WR, WI, VS, &ldvs,
	  work, &lwork, bwork, info);
    lwork = (int) work[0];
    work = (double *) realloc (work, lwork * sizeof(double));
    dgees_(&jobv, &sort, select, &n, A, &lda, sdim, WR, WI, VS, &ldvs,
	  work, &lwork, bwork, info);
    free(bwork);
    free(work);
}
//void dgeesx(char, char, int (*)(double *, double *), char, int,
//	    double *, int, int *, double *, double *, double *,
//	    int, double *, double *, int *) {
//}
void dgeev(char jobvl, char jobvr, int n, double *A, int lda,
	   double *WR, double *WI, double *VL, int ldvl,
	   double *VR, int ldvr, int *info) {
    int lwork = -1;
    double *work = (double *) malloc (sizeof(double));
    // query workspace size
    dgeev_(&jobvl, &jobvr, &n, A, &lda, WR, WI, VL, &ldvl,
	   VR, &ldvr, work, &lwork, info);
    lwork = (int) work[0];
    work = (double *) realloc (work, lwork * sizeof(double));
    dgeev_(&jobvl, &jobvr, &n, A, &lda, WR, WI, VL, &ldvl,
	   VR, &ldvr, work, &lwork, info);
    free(work);
}
//void dgeevx(char, char, char, char, int, double *, int, double *,
//	    double *, double *, int, double *, int, int *, int *,
//	    double *, double *, double *, double *, int *) {
//}
void zgeev(char jobvl, char jobvr, int n, Complex *A, int lda,
	   Complex *W, Complex *VL, int ldvl, Complex *VR, int ldvr,
	   int *info) {
    int lwork = -1;
    Complex *work = (Complex *) malloc (sizeof(Complex));
    double *rwork = (double *) malloc (2*n*sizeof(double));
    // query workspace size
    zgeev_(&jobvl, &jobvr, &n, A, &lda, W, VL, &ldvl,
	   VR, &ldvr, work, &lwork, rwork, info);
    lwork = (int) *((double *) work);
    work = (Complex *) realloc (work, lwork * sizeof(Complex));
    zgeev_(&jobvl, &jobvr, &n, A, &lda, W, VL, &ldvl,
	   VR, &ldvr, work, &lwork, rwork, info);
    free(rwork);
    free(work);
}

void dgesvd(char jobu, char jobvt, int m, int n, double *A, int lda,
	    double *S, double *U, int ldu, double *VT, int ldvt,
	    int *info) {
    int lwork = -1;
    double *work = (double *) malloc (sizeof(double));
    // query workspace size
    dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt,
	    work, &lwork, info);
    lwork = (int) *((double *) work);
    work = (double *) realloc (work, lwork * sizeof(double));
    dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt,
	    work, &lwork, info);
    free(work);
}

void zgesvd(char jobu, char jobvt, int m, int n, Complex *A, int lda,
	    double *S, Complex *U, int ldu, Complex *VT, int ldvt,
	    int *info) {
    int lwork = -1;
    Complex *work = (Complex *) malloc (sizeof(Complex));
    double *rwork = (double *) malloc (5*n*sizeof(double));
    // query workspace size
    zgesvd_(&jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt,
	    work, &lwork, rwork, info);
    lwork = (int) *((double *) work);
    work = (Complex *) realloc (work, lwork * sizeof(Complex));
    zgesvd_(&jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt,
	    work, &lwork, rwork, info);
    free(rwork);
    free(work);
}

void dsyev(char jobz, char uplo, int n, double *A, int lda,
	   double *W, int *info) {
    int lwork = -1;
    double *work = (double *) malloc (sizeof(double));
    // query workspace size
    dsyev_(&jobz, &uplo, &n, A, &lda, W, work, &lwork, info);
    lwork = (int) work[0];
    work = (double *) realloc (work, lwork * sizeof(double));
    dsyev_(&jobz, &uplo, &n, A, &lda, W, work, &lwork, info);
    free(work);
}

void zheev(char jobz, char uplo, int n, Complex *A, int lda,
	   double *W, int *info) {
    int lwork = -1;
    Complex *work = (Complex *) malloc (sizeof(Complex));
    double *rwork = (double *) malloc ((3*n+1)*sizeof(double));
    // query workspace size
    zheev_(&jobz, &uplo, &n, A, &lda, W, work, &lwork, rwork, info);
    lwork = (int) *((double *) work);
    work = (Complex *) realloc (work, lwork * sizeof(Complex));
    zheev_(&jobz, &uplo, &n, A, &lda, W, work, &lwork, rwork, info);
    free(rwork);
    free(work);
}
