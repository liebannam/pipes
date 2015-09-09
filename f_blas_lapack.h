#ifndef __CLAPACK_H
#define __CLAPACK_H

#ifndef JON_COMPLEX
#define JON_COMPLEX
typedef struct single_cx_struct { float r, i; } single_cx;
typedef struct double_cx_struct { double r, i; } double_cx;
#endif

typedef int (*L_fp)();

float 
sdot_(int* N, 
      float* X, int* incX, 
      float* Y, int* incY);

double
ddot_(int* N, 
      double* X, int* incX, 
      double* Y, int* incY);

void 
cdotu_(single_cx* retval,
       int* N, 
       single_cx* X, int* incX, 
       single_cx* Y, int* incY);

void
cdotc_(single_cx* retval,
       int* N, 
       single_cx* X, int* incX, 
       single_cx* Y, int* incY);

void
zdotu_(double_cx* retval,
       int* N, 
       double_cx* X, int* incX, 
       double_cx* Y, int* incY);

void
zdotc_(double_cx* retval,
       int* N, 
       double_cx* X, int* incX, 
       double_cx* Y, int* incY);

float 
snrm2_(int* N, 
       float* X, int* incX);

float
sasum_(int* N, 
       float* X, int* incX);

double
dnrm2_(int* N, 
       double* X, int* incX);

double
dasum_(int* N, 
       double* X, int* incX);

float 
scnrm2_(int* N, 
        single_cx* X, int* incX);

float
scasum_(int* N, 
        single_cx* X, int* incX);

double 
dznrm2_(int* N, 
        double_cx* X, int* incX);

double
dzasum_(int* N, 
        double_cx* X, int* incX);

int
isamax_(int* N,
        float* X, int* incX);

int
idamax_(int* N,
        double* X, int* incX);

int
icamax_(int* N,
        single_cx* X, int* incX);

int
izamax_(int* N,
        double_cx* X, int* incX);

int
sswap_(int* N,
       float* X, int* incX,
       float* Y, int* incY);

int
scopy_(int* N,
       float* X, int* incX,
       float* Y, int* incY);

int
saxpy_(int* N,
       float* alpha,
       float* X, int* incX,
       float* Y, int* incY);

int
dswap_(int* N,
       double* X, int* incX,
       double* Y, int* incY);

int
dcopy_(int* N,
       double* X, int* incX,
       double* Y, int* incY);

int
daxpy_(int* N,
       double* alpha,
       double* X, int* incX,
       double* Y, int* incY);

int
cswap_(int* N,
       single_cx* X, int* incX,
       single_cx* Y, int* incY);

int
ccopy_(int* N,
       single_cx* X, int* incX,
       single_cx* Y, int* incY);

int
caxpy_(int* N,
      single_cx* alpha,
      single_cx* X, int* incX,
      single_cx* Y, int* incY);

int
zswap_(int* N,
       double_cx* X, int* incX,
       double_cx* Y, int* incY);

int
zcopy_(int* N,
       double_cx* X, int* incX,
       double_cx* Y, int* incY);

int
zaxpy_(int* N,
       double_cx* alpha,
       double_cx* X, int* incX,
       double_cx* Y, int* incY);

int
srotg_(float* a, float* b, float* c, float* s);

int
srot_(int* N,
      float* X, int* incX,
      float* Y, int* incY,
      float* c, float* s);

int
drotg_(double* a, double* b, double* c, double* s);

int
drot_(int* N,
      double* X, int* incX,
      double* Y, int* incY,
      double* c, double* s);

int
sscal_(int* N,
       float* alpha,
       float* X, int* incX);

int
dscal_(int* N,
       double* alpha,
       double* X, int* incX);

int
cscal_(int* N,
       single_cx* alpha,
       single_cx* X, int* incX);

int
zscal_(int* N,
       double_cx* alpha,
       double_cx* X, int* incX);

int
csscal_(int* N,
        float* alpha,
        single_cx* X, int* incX);

int
zdscal_(int* N,
        double* alpha,
        double_cx* X, int* incX);

int
sgemv_(char* trans, int* M, int* N,
       float* alpha,
       float* A, int* lda,
       float* X, int* incX,
       float* beta,
       float* Y, int* incY);

int
sgbmv_(char *trans, int *M, int *N, int *KL, int *KU, 
       float *alpha, 
       float *A, int *lda, 
       float *X, int *incX, 
       float *beta, 
       float *Y, int *incY);

int 
strmv_(char* uplo, char *trans, char* diag, int *N,  
       float *A, int *lda, 
       float *X, int *incX);

int
stbmv_(char* uplo, char* trans, char* diag, int* N, int* K,
       float* A, int* lda,
       float* X, int* incX);

int
stpmv_(char* uplo, char* trans, char* diag, int* N, 
       float* Ap, 
       float* X, int* incX);

int
strsv_(char* uplo, char* trans, char* diag, int* N,
       float* A, int* lda,
       float* X, int* incX);

int
stbsv_(char* uplo, char* trans, char* diag, int* N, int* K,
       float* A, int* lda, 
       float* X, int* incX);

int
stpsv_(char* uplo, char* trans, char* diag, int* N, 
       float* Ap, 
       float* X, int* incX);

int
dgemv_(char* trans, int* M, int* N,
       double* alpha,
       double* A, int* lda,
       double* X, int* incX,
       double* beta,
       double* Y, int* incY);

int 
dgbmv_(char *trans, int *M, int *N, int *KL, int *KU, 
       double *alpha, 
       double *A, int *lda, 
       double *X, int *incX, 
       double *beta, 
       double *Y, int *incY);

int 
dtrmv_(char* uplo, char *trans, char* diag, int *N,  
       double *A, int *lda, 
       double *X, int *incX);

int
dtbmv_(char* uplo, char* trans, char* diag, int* N, int* K,
       double* A, int* lda,
       double* X, int* incX);

int
dtpmv_(char* uplo, char* trans, char* diag, int* N, 
       double* Ap, 
       double* X, int* incX);

int
dtrsv_(char* uplo, char* trans, char* diag, int* N,
       double* A, int* lda,
       double* X, int* incX);

int
dtbsv_(char* uplo, char* trans, char* diag, int* N, int* K,
       double* A, int* lda, 
       double* X, int* incX);

int
dtpsv_(char* uplo, char* trans, char* diag, int* N, 
       double* Ap, 
       double* X, int* incX);

int
cgemv_(char* trans, int* M, int* N,
       single_cx* alpha,
       single_cx* A, int* lda,
       single_cx* X, int* incX,
       single_cx* beta,
       single_cx* Y, int* incY);

int 
cgbmv_(char *trans, int *M, int *N, int *KL, int *KU, 
       single_cx *alpha, 
       single_cx *A, int *lda, 
       single_cx *X, int *incX, 
       single_cx *beta, 
       single_cx *Y, int *incY);

int 
ctrmv_(char* uplo, char *trans, char* diag, int *N,  
       single_cx *A, int *lda, 
       single_cx *X, int *incX);

int
ctbmv_(char* uplo, char* trans, char* diag, int* N, int* K,
       single_cx* A, int* lda,
       single_cx* X, int* incX);

int
ctpmv_(char* uplo, char* trans, char* diag, int* N, 
       single_cx* Ap, 
       single_cx* X, int* incX);

int
ctrsv_(char* uplo, char* trans, char* diag, int* N,
       single_cx* A, int* lda,
       single_cx* X, int* incX);

int
ctbsv_(char* uplo, char* trans, char* diag, int* N, int* K,
       single_cx* A, int* lda, 
       single_cx* X, int* incX);

int
ctpsv_(char* uplo, char* trans, char* diag, int* N, 
       single_cx* Ap, 
       single_cx* X, int* incX);

int
zgemv_(char* trans, int* M, int* N,
       double_cx* alpha,
       double_cx* A, int* lda,
       double_cx* X, int* incX,
       double_cx* beta,
       double_cx* Y, int* incY);

int 
zgbmv_(char *trans, int *M, int *N, int *KL, int *KU, 
       double_cx *alpha, 
       double_cx *A, int *lda, 
       double_cx *X, int *incX, 
       double_cx *beta, 
       double_cx *Y, int *incY);

int 
ztrmv_(char* uplo, char *trans, char* diag, int *N,  
       double_cx *A, int *lda, 
       double_cx *X, int *incX);

int
ztbmv_(char* uplo, char* trans, char* diag, int* N, int* K,
       double_cx* A, int* lda,
       double_cx* X, int* incX);

 void  
ztpmv_(char* uplo, char* trans, char* diag, int* N, 
      double_cx* Ap, 
      double_cx* X, int* incX);

int
ztrsv_(char* uplo, char* trans, char* diag, int* N,
       double_cx* A, int* lda,
       double_cx* X, int* incX);

int
ztbsv_(char* uplo, char* trans, char* diag, int* N, int* K,
       double_cx* A, int* lda, 
       double_cx* X, int* incX);

int
ztpsv_(char* uplo, char* trans, char* diag, int* N, 
       double_cx* Ap, 
       double_cx* X, int* incX);

int
ssymv_(char* uplo, int* N,
       float* alpha,
       float* A, int* lda,
       float* X, int* incX,
       float* beta,
       float* Y, int* incY);

int 
ssbmv_(char* uplo, int* N, int* K,
       float* alpha,
       float* A, int* lda,
       float* X, int* incX,
       float* beta,
       float* Y, int* incY);

int
sspmv_(char* uplo, int* N,
       float* alpha,
       float* Ap,
       float* X, int* incX,
       float* beta,
       float* Y, int* incY);

int
sger_(int* M, int* N,
      float* alpha,
      float* X, int* incX,
      float* Y, int* incY,
      float* A, int* lda);

int
ssyr_(char* uplo, int* N,
      float* alpha,
      float* X, int* incX,
      float* A, int* lda);

int
sspr_(char* uplo, int* N,
      float* alpha,
      float* X, int* incX,
      float* Ap);

int
ssyr2_(char* uplo, int* N,
       float* alpha,
       float* X, int* incX,
       float* Y, int* incY,
       float* A, int* lda);

int
sspr2_(char* uplo, int* N,
       float* alpha, 
       float* X, int* incX,
       float* Y, int* incY,
       float* A);

int
dsymv_(char* uplo, int* N,
       double* alpha,
       double* A, int* lda,
       double* X, int* incX,
       double* beta,
       double* Y, int* incY);

int 
dsbmv_(char* uplo, int* N, int* K,
       double* alpha,
       double* A, int* lda,
       double* X, int* incX,
       double* beta,
       double* Y, int* incY);

int
dspmv_(char* uplo, int* N,
       double* alpha,
       double* Ap,
       double* X, int* incX,
       double* beta,
       double* Y, int* incY);

int
dger_(int* M, int* N,
      double* alpha,
      double* X, int* incX,
      double* Y, int* incY,
      double* A, int* lda);

int
dsyr_(char* uplo, int* N,
      double* alpha,
      double* X, int* incX,
      double* A, int* lda);

int
dspr_(char* uplo, int* N,
      double* alpha,
      double* X, int* incX,
      double* Ap);

int
dsyr2_(char* uplo, int* N,
       double* alpha,
       double* X, int* incX,
       double* Y, int* incY,
       double* A, int* lda);

int
dspr2_(char* uplo, int* N,
       double* alpha, 
       double* X, int* incX,
       double* Y, int* incY,
       double* A);

int
chemv_(char* uplo, int* N,
       single_cx* alpha,
       single_cx* A, int* lda,
       single_cx* X, int* incX,
       single_cx* beta,
       single_cx* Y, int* incY);

int
chbmv_(char* uplo, int* N, int* K,
       single_cx* alpha,
       single_cx* A, int* lda,
       single_cx* X, int* incX,
       single_cx* beta,
       single_cx* Y, int* incY);

int
chpmv_(char* uplo, int* N, 
       single_cx* alpha,
       single_cx* Ap, 
       single_cx* X, int* incX,
       single_cx* beta,
       single_cx* Y, int* incY);

int
cgeru_(int* M, int* N,
       single_cx* alpha,
       single_cx* X, int* incX,
       single_cx* Y, int* incY,
       single_cx* A, int* lda);

int
cgerc_(int* M, int* N,
       single_cx* alpha,
       single_cx* X, int* incX,
       single_cx* Y, int* incY,
       single_cx* A, int* lda);

int
cher_(char* uplo, int* N,
      float* alpha,
      single_cx* X, int* incX,
      single_cx* A, int* lda);

int
chpr_(char* uplo, int* N,
      float* alpha,
      single_cx* X, int* incX,
      single_cx* Ap);

int
cher2_(char* uplo, int* N,
       single_cx* alpha,
       single_cx* X, int* incX,
       single_cx* Y, int* incY,
       single_cx* A, int* lda);

int
chpr2_(char* uplo, int* N,
       single_cx* alpha,
       single_cx* X, int* incX,
       single_cx* Y, int* incY,
       single_cx* Ap);

int
zhemv_(char* uplo, int* N,
       double_cx* alpha,
       double_cx* A, int* lda,
       double_cx* X, int* incX,
       double_cx* beta,
       double_cx* Y, int* incY);

int
zhbmv_(char* uplo, int* N, int* K,
       double_cx* alpha,
       double_cx* A, int* lda,
       double_cx* X, int* incX,
       double_cx* beta,
       double_cx* Y, int* incY);

int
zhpmv_(char* uplo, int* N, 
       double_cx* alpha,
       double_cx* Ap, 
       double_cx* X, int* incX,
       double_cx* beta,
       double_cx* Y, int* incY);

int
zgeru_(int* M, int* N,
       double_cx* alpha,
       double_cx* X, int* incX,
       double_cx* Y, int* incY,
       double_cx* A, int* lda);

int
zgerc_(int* M, int* N,
       double_cx* alpha,
       double_cx* X, int* incX,
       double_cx* Y, int* incY,
       double_cx* A, int* lda);

int
zher_(char* uplo, int* N,
      double* alpha,
      double_cx* X, int* incX,
      double_cx* A, int* lda);

int
zhpr_(char* uplo, int* N,
      double* alpha,
      double_cx* X, int* incX,
      double_cx* Ap);

int
zher2_(char* uplo, int* N,
       double_cx* alpha,
       double_cx* X, int* incX,
       double_cx* Y, int* incY,
       double_cx* A, int* lda);

int
zhpr2_(char* uplo, int* N,
       double_cx* alpha,
       double_cx* X, int* incX,
       double_cx* Y, int* incY,
       double_cx* Ap);

int
sgemm_(char* transA, char* transB, int* M, int* N, int* K,
       float* alpha,
       float* A, int* lda,
       float* B, int* ldb,
       float* beta,
       float* C, int* ldc);

int
ssymm_(char* side, char* uplo, int* M, int* N,
       float* alpha,
       float* A, int* lda,
       float* B, int* ldb,
       float* beta,
       float* C, int* ldc);

int
ssyrk_(char* uplo, char* trans, int* N, int* K,
       float* alpha,
       float* A, int* lda,
       float* beta,
       float* C, int* ldc);

int
ssyr2k_(char* uplo, char* trans, int* N, int* K,
        float* alpha,
        float* A, int* lda,
        float* B, int* ldb,
        float* beta,
        float* C, int* ldc);

int
strmm_(char* side, char* uplo, char* trans, char* diag, 
       int* M, int* N,
       float* alpha,
       float* A, int* lda,
       float* B, int* ldb);

int 
strsm_(char* side, char* uplo, char* trans, char* diag,
       int* M, int* N,
       float* alpha,
       float* A, int* lda,
       float* B, int* ldb);

int
dgemm_(char* transA, char* transB, int* M, int* N, int* K,
       double* alpha,
       double* A, int* lda,
       double* B, int* ldb,
       double* beta,
       double* C, int* ldc);

int
dsymm_(char* side, char* uplo, int* M, int* N,
       double* alpha,
       double* A, int* lda,
       double* B, int* ldb,
       double* beta,
       double* C, int* ldc);

int
dsyrk_(char* uplo, char* trans, int* N, int* K,
       double* alpha,
       double* A, int* lda,
       double* beta,
       double* C, int* ldc);

int
dsyr2k_(char* uplo, char* trans, int* N, int* K,
        double* alpha,
        double* A, int* lda,
        double* B, int* ldb,
        double* beta,
        double* C, int* ldc);

int
dtrmm_(char* side, char* uplo, char* trans, char* diag, 
       int* M, int* N,
       double* alpha,
       double* A, int* lda,
       double* B, int* ldb);

int 
dtrsm_(char* side, char* uplo, char* trans, char* diag,
       int* M, int* N,
       double* alpha,
       double* A, int* lda,
       double* B, int* ldb);

int
cgemm_(char* transA, char* transB, int* M, int* N, int* K,
       single_cx* alpha,
       single_cx* A, int* lda,
       single_cx* B, int* ldb,
       single_cx* beta,
       single_cx* C, int* ldc);

int
csymm_(char* side, char* uplo, int* M, int* N,
       single_cx* alpha,
       single_cx* A, int* lda,
       single_cx* B, int* ldb,
       single_cx* beta,
       single_cx* C, int* ldc);

int
csyrk_(char* uplo, char* trans, int* N, int* K,
       single_cx* alpha,
       single_cx* A, int* lda,
       single_cx* beta,
       single_cx* C, int* ldc);

int
csyr2k_(char* uplo, char* trans, int* N, int* K,
        single_cx* alpha,
        single_cx* A, int* lda,
        single_cx* B, int* ldb,
        single_cx* beta,
        single_cx* C, int* ldc);

int
ctrmm_(char* side, char* uplo, char* trans, char* diag, 
       int* M, int* N,
       single_cx* alpha,
       single_cx* A, int* lda,
       single_cx* B, int* ldb);

int 
ctrsm_(char* side, char* uplo, char* trans, char* diag,
       int* M, int* N,
       single_cx* alpha,
       single_cx* A, int* lda,
       single_cx* B, int* ldb);

int
zgemm_(char* transA, char* transB, int* M, int* N, int* K,
       double_cx* alpha,
       double_cx* A, int* lda,
       double_cx* B, int* ldb,
       double_cx* beta,
       double_cx* C, int* ldc);

int
zsymm_(char* side, char* uplo, int* M, int* N,
       double_cx* alpha,
       double_cx* A, int* lda,
       double_cx* B, int* ldb,
       double_cx* beta,
       double_cx* C, int* ldc);

int
zsyrk_(char* uplo, char* trans, int* N, int* K,
       double_cx* alpha,
       double_cx* A, int* lda,
       double_cx* beta,
       double_cx* C, int* ldc);

int
zsyr2k_(char* uplo, char* trans, int* N, int* K,
        double_cx* alpha,
        double_cx* A, int* lda,
        double_cx* B, int* ldb,
        double_cx* beta,
        double_cx* C, int* ldc);

int
ztrmm_(char* side, char* uplo, char* trans, char* diag, 
       int* M, int* N,
       double_cx* alpha,
       double_cx* A, int* lda,
       double_cx* B, int* ldb);

int 
ztrsm_(char* side, char* uplo, char* trans, char* diag,
       int* M, int* N,
       double_cx* alpha,
       double_cx* A, int* lda,
       double_cx* B, int* ldb);

int
chemm_(char* side, char* uplo, int* M, int* N,
       single_cx* alpha,
       single_cx* A, int* lda,
       single_cx* B, int* ldb,
       single_cx* beta,
       single_cx* C, int* ldc);

int
cherk_(char* uplo, char* trans, int* N, int* K,
       float* alpha,
       single_cx* A, int* lda,
       float* beta,
       single_cx* C, int* ldc);

int
cher2k_(char* uplo, char* trans, int* N, int* K,
        single_cx* alpha,
        single_cx* A, int* lda,
        single_cx* B, int* ldb,
        float* beta,
        single_cx* C, int* ldc);

int
zhemm_(char* side, char* uplo, int* M, int* N,
       double_cx* alpha,
       double_cx* A, int* lda,
       double_cx* B, int* ldb,
       double_cx* beta,
       double_cx* C, int* ldc);

int
zherk_(char* uplo, char* trans, int* N, int* K,
       double* alpha,
       double_cx* A, int* lda,
       double* beta,
       double_cx* C, int* ldc);

int
zher2k_(char* uplo, char* trans, int* N, int* K,
        double_cx* alpha,
        double_cx* A, int* lda,
        double_cx* B, int* ldb,
        double* beta,
        double_cx* C, int* ldc);
 
int cbdsqr_(char *uplo, int *n, int *ncvt, int *
	nru, int *ncc, float *d__, float *e, single_cx *vt, int *ldvt, 
	single_cx *u, int *ldu, single_cx *c__, int *ldc, float *rwork, 
	int *info);
 
int cgbbrd_(char *vect, int *m, int *n, int *ncc,
	 int *kl, int *ku, single_cx *ab, int *ldab, float *d__, 
	float *e, single_cx *q, int *ldq, single_cx *pt, int *ldpt, 
	single_cx *c__, int *ldc, single_cx *work, float *rwork, int *info);
 
int cgbcon_(char *norm, int *n, int *kl, int *ku,
	 single_cx *ab, int *ldab, int *ipiv, float *anorm, float *rcond, 
	single_cx *work, float *rwork, int *info);
 
int cgbequ_(int *m, int *n, int *kl, int *ku,
	 single_cx *ab, int *ldab, float *r__, float *c__, float *rowcnd, float 
	*colcnd, float *amax, int *info);
 
int cgbrfs_(char *trans, int *n, int *kl, int *
	ku, int *nrhs, single_cx *ab, int *ldab, single_cx *afb, int *
	ldafb, int *ipiv, single_cx *b, int *ldb, single_cx *x, int *
	ldx, float *ferr, float *berr, single_cx *work, float *rwork, int *
	info);
 
int cgbsv_(int *n, int *kl, int *ku, int *
	nrhs, single_cx *ab, int *ldab, int *ipiv, single_cx *b, int *
	ldb, int *info);
 
int cgbsvx_(char *fact, char *trans, int *n, int *kl,
	 int *ku, int *nrhs, single_cx *ab, int *ldab, single_cx *afb,
	 int *ldafb, int *ipiv, char *equed, float *r__, float *c__, 
	single_cx *b, int *ldb, single_cx *x, int *ldx, float *rcond, float 
	*ferr, float *berr, single_cx *work, float *rwork, int *info);
 
int cgbtf2_(int *m, int *n, int *kl, int *ku,
	 single_cx *ab, int *ldab, int *ipiv, int *info);
 
int cgbtrf_(int *m, int *n, int *kl, int *ku,
	 single_cx *ab, int *ldab, int *ipiv, int *info);
 
int cgbtrs_(char *trans, int *n, int *kl, int *
	ku, int *nrhs, single_cx *ab, int *ldab, int *ipiv, single_cx 
	*b, int *ldb, int *info);
 
int cgebak_(char *job, char *side, int *n, int *ilo, 
	int *ihi, float *scale, int *m, single_cx *v, int *ldv, 
	int *info);
 
int cgebal_(char *job, int *n, single_cx *a, int *lda, 
	int *ilo, int *ihi, float *scale, int *info);
 
int cgebd2_(int *m, int *n, single_cx *a, int *lda,
	 float *d__, float *e, single_cx *tauq, single_cx *taup, single_cx *work, 
	int *info);
 
int cgebrd_(int *m, int *n, single_cx *a, int *lda,
	 float *d__, float *e, single_cx *tauq, single_cx *taup, single_cx *work, 
	int *lwork, int *info);
 
int cgecon_(char *norm, int *n, single_cx *a, int *lda,
	 float *anorm, float *rcond, single_cx *work, float *rwork, int *info);
 
int cgeequ_(int *m, int *n, single_cx *a, int *lda,
	 float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, 
	int *info);
 
int cgees_(char *jobvs, char *sort, L_fp select, int *n, 
	single_cx *a, int *lda, int *sdim, single_cx *w, single_cx *vs, 
	int *ldvs, single_cx *work, int *lwork, float *rwork, int *
	bwork, int *info);
 
int cgeesx_(char *jobvs, char *sort, L_fp select, char *
	sense, int *n, single_cx *a, int *lda, int *sdim, single_cx *
	w, single_cx *vs, int *ldvs, float *rconde, float *rcondv, single_cx *
	work, int *lwork, float *rwork, int *bwork, int *info);
 
int cgeev_(char *jobvl, char *jobvr, int *n, single_cx *a, 
	int *lda, single_cx *w, single_cx *vl, int *ldvl, single_cx *vr, 
	int *ldvr, single_cx *work, int *lwork, float *rwork, int *
	info);
 
int cgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, int *n, single_cx *a, int *lda, single_cx *w, single_cx *vl, 
	int *ldvl, single_cx *vr, int *ldvr, int *ilo, int *ihi,
	 float *scale, float *abnrm, float *rconde, float *rcondv, single_cx *work, 
	int *lwork, float *rwork, int *info);
 
int cgegs_(char *jobvsl, char *jobvsr, int *n, single_cx *
	a, int *lda, single_cx *b, int *ldb, single_cx *alpha, single_cx *
	beta, single_cx *vsl, int *ldvsl, single_cx *vsr, int *ldvsr, 
	single_cx *work, int *lwork, float *rwork, int *info);
 
int cgegv_(char *jobvl, char *jobvr, int *n, single_cx *a, 
	int *lda, single_cx *b, int *ldb, single_cx *alpha, single_cx *beta,
	 single_cx *vl, int *ldvl, single_cx *vr, int *ldvr, single_cx *
	work, int *lwork, float *rwork, int *info);
 
int cgehd2_(int *n, int *ilo, int *ihi, single_cx *
	a, int *lda, single_cx *tau, single_cx *work, int *info);
 
int cgehrd_(int *n, int *ilo, int *ihi, single_cx *
	a, int *lda, single_cx *tau, single_cx *work, int *lwork, int 
	*info);
 
int cgelq2_(int *m, int *n, single_cx *a, int *lda,
	 single_cx *tau, single_cx *work, int *info);
 
int cgelqf_(int *m, int *n, single_cx *a, int *lda,
	 single_cx *tau, single_cx *work, int *lwork, int *info);
 
int cgels_(char *trans, int *m, int *n, int *
	nrhs, single_cx *a, int *lda, single_cx *b, int *ldb, single_cx *
	work, int *lwork, int *info);
 
int cgelsx_(int *m, int *n, int *nrhs, single_cx *
	a, int *lda, single_cx *b, int *ldb, int *jpvt, float *rcond,
	 int *rank, single_cx *work, float *rwork, int *info);
 
int cgelsy_(int *m, int *n, int *nrhs, single_cx *
	a, int *lda, single_cx *b, int *ldb, int *jpvt, float *rcond,
	 int *rank, single_cx *work, int *lwork, float *rwork, int *
	info);
 
int cgeql2_(int *m, int *n, single_cx *a, int *lda,
	 single_cx *tau, single_cx *work, int *info);
 
int cgeqlf_(int *m, int *n, single_cx *a, int *lda,
	 single_cx *tau, single_cx *work, int *lwork, int *info);
 
int cgeqp3_(int *m, int *n, single_cx *a, int *lda,
	 int *jpvt, single_cx *tau, single_cx *work, int *lwork, float *
	rwork, int *info);
 
int cgeqpf_(int *m, int *n, single_cx *a, int *lda,
	 int *jpvt, single_cx *tau, single_cx *work, float *rwork, int *
	info);
 
int cgeqr2_(int *m, int *n, single_cx *a, int *lda,
	 single_cx *tau, single_cx *work, int *info);
 
int cgeqrf_(int *m, int *n, single_cx *a, int *lda,
	 single_cx *tau, single_cx *work, int *lwork, int *info);
 
int cgerfs_(char *trans, int *n, int *nrhs, single_cx *
	a, int *lda, single_cx *af, int *ldaf, int *ipiv, single_cx *
	b, int *ldb, single_cx *x, int *ldx, float *ferr, float *berr, 
	single_cx *work, float *rwork, int *info);
 
int cgerq2_(int *m, int *n, single_cx *a, int *lda,
	 single_cx *tau, single_cx *work, int *info);
 
int cgerqf_(int *m, int *n, single_cx *a, int *lda,
	 single_cx *tau, single_cx *work, int *lwork, int *info);
 
int cgesc2_(int *n, single_cx *a, int *lda, single_cx *
	rhs, int *ipiv, int *jpiv, float *scale);
 
int cgesv_(int *n, int *nrhs, single_cx *a, int *
	lda, int *ipiv, single_cx *b, int *ldb, int *info);
 
int cgesvx_(char *fact, char *trans, int *n, int *
	nrhs, single_cx *a, int *lda, single_cx *af, int *ldaf, int *
	ipiv, char *equed, float *r__, float *c__, single_cx *b, int *ldb, 
	single_cx *x, int *ldx, float *rcond, float *ferr, float *berr, 
	single_cx *work, float *rwork, int *info);
 
int cgetc2_(int *n, single_cx *a, int *lda, int *
	ipiv, int *jpiv, int *info);
 
int cgetf2_(int *m, int *n, single_cx *a, int *lda,
	 int *ipiv, int *info);
 
int cgetrf_(int *m, int *n, single_cx *a, int *lda,
	 int *ipiv, int *info);
 
int cgetri_(int *n, single_cx *a, int *lda, int *
	ipiv, single_cx *work, int *lwork, int *info);
 
int cgetrs_(char *trans, int *n, int *nrhs, single_cx *
	a, int *lda, int *ipiv, single_cx *b, int *ldb, int *
	info);
 
int cggbak_(char *job, char *side, int *n, int *ilo, 
	int *ihi, float *lscale, float *rscale, int *m, single_cx *v, 
	int *ldv, int *info);
 
int cggbal_(char *job, int *n, single_cx *a, int *lda, 
	single_cx *b, int *ldb, int *ilo, int *ihi, float *lscale, 
	float *rscale, float *work, int *info);
 
int cgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, int *n, single_cx *a, int *lda, single_cx *b, int *
	ldb, int *sdim, single_cx *alpha, single_cx *beta, single_cx *vsl, 
	int *ldvsl, single_cx *vsr, int *ldvsr, single_cx *work, int *
	lwork, float *rwork, int *bwork, int *info);
 
int cggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, char *sense, int *n, single_cx *a, int *lda, single_cx *b,
	 int *ldb, int *sdim, single_cx *alpha, single_cx *beta, single_cx *
	vsl, int *ldvsl, single_cx *vsr, int *ldvsr, float *rconde, float 
	*rcondv, single_cx *work, int *lwork, float *rwork, int *iwork, 
	int *liwork, int *bwork, int *info);
 
int cggev_(char *jobvl, char *jobvr, int *n, single_cx *a, 
	int *lda, single_cx *b, int *ldb, single_cx *alpha, single_cx *beta,
	 single_cx *vl, int *ldvl, single_cx *vr, int *ldvr, single_cx *
	work, int *lwork, float *rwork, int *info);
 
int cggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, int *n, single_cx *a, int *lda, single_cx *b, int *ldb,
	 single_cx *alpha, single_cx *beta, single_cx *vl, int *ldvl, single_cx *
	vr, int *ldvr, int *ilo, int *ihi, float *lscale, float *
	rscale, float *abnrm, float *bbnrm, float *rconde, float *rcondv, single_cx 
	*work, int *lwork, float *rwork, int *iwork, int *bwork, 
	int *info);
 
int cggglm_(int *n, int *m, int *p, single_cx *a, 
	int *lda, single_cx *b, int *ldb, single_cx *d__, single_cx *x, 
	single_cx *y, single_cx *work, int *lwork, int *info);
 
int cgghrd_(char *compq, char *compz, int *n, int *
	ilo, int *ihi, single_cx *a, int *lda, single_cx *b, int *ldb,
	 single_cx *q, int *ldq, single_cx *z__, int *ldz, int *info);
 
int cgglse_(int *m, int *n, int *p, single_cx *a, 
	int *lda, single_cx *b, int *ldb, single_cx *c__, single_cx *d__, 
	single_cx *x, single_cx *work, int *lwork, int *info);
 
int cggqrf_(int *n, int *m, int *p, single_cx *a, 
	int *lda, single_cx *taua, single_cx *b, int *ldb, single_cx *taub, 
	single_cx *work, int *lwork, int *info);
 
int cggrqf_(int *m, int *p, int *n, single_cx *a, 
	int *lda, single_cx *taua, single_cx *b, int *ldb, single_cx *taub, 
	single_cx *work, int *lwork, int *info);
 
int cggsvd_(char *jobu, char *jobv, char *jobq, int *m, 
	int *n, int *p, int *k, int *l, single_cx *a, int *
	lda, single_cx *b, int *ldb, float *alpha, float *beta, single_cx *u, 
	int *ldu, single_cx *v, int *ldv, single_cx *q, int *ldq, 
	single_cx *work, float *rwork, int *iwork, int *info);
 
int cggsvp_(char *jobu, char *jobv, char *jobq, int *m, 
	int *p, int *n, single_cx *a, int *lda, single_cx *b, int 
	*ldb, float *tola, float *tolb, int *k, int *l, single_cx *u, 
	int *ldu, single_cx *v, int *ldv, single_cx *q, int *ldq, 
	int *iwork, float *rwork, single_cx *tau, single_cx *work, int *
	info);
 
int cgtcon_(char *norm, int *n, single_cx *dl, single_cx *
	d__, single_cx *du, single_cx *du2, int *ipiv, float *anorm, float *
	rcond, single_cx *work, int *info);
 
int cgtrfs_(char *trans, int *n, int *nrhs, single_cx *
	dl, single_cx *d__, single_cx *du, single_cx *dlf, single_cx *df, single_cx *
	duf, single_cx *du2, int *ipiv, single_cx *b, int *ldb, single_cx *
	x, int *ldx, float *ferr, float *berr, single_cx *work, float *rwork, 
	int *info);
 
int cgtsv_(int *n, int *nrhs, single_cx *dl, single_cx *
	d__, single_cx *du, single_cx *b, int *ldb, int *info);
 
int cgtsvx_(char *fact, char *trans, int *n, int *
	nrhs, single_cx *dl, single_cx *d__, single_cx *du, single_cx *dlf, single_cx *
	df, single_cx *duf, single_cx *du2, int *ipiv, single_cx *b, int *
	ldb, single_cx *x, int *ldx, float *rcond, float *ferr, float *berr, 
	single_cx *work, float *rwork, int *info);
 
int cgttrf_(int *n, single_cx *dl, single_cx *d__, single_cx *
	du, single_cx *du2, int *ipiv, int *info);
 
int cgttrs_(char *trans, int *n, int *nrhs, single_cx *
	dl, single_cx *d__, single_cx *du, single_cx *du2, int *ipiv, single_cx *
	b, int *ldb, int *info);
 
int cgtts2_(int *itrans, int *n, int *nrhs, 
	single_cx *dl, single_cx *d__, single_cx *du, single_cx *du2, int *ipiv, 
	single_cx *b, int *ldb);
 
int chbev_(char *jobz, char *uplo, int *n, int *kd, 
	single_cx *ab, int *ldab, float *w, single_cx *z__, int *ldz, 
	single_cx *work, float *rwork, int *info);
 
int chbevd_(char *jobz, char *uplo, int *n, int *kd, 
	single_cx *ab, int *ldab, float *w, single_cx *z__, int *ldz, 
	single_cx *work, int *lwork, float *rwork, int *lrwork, int *
	iwork, int *liwork, int *info);
 
int chbevx_(char *jobz, char *range, char *uplo, int *n, 
	int *kd, single_cx *ab, int *ldab, single_cx *q, int *ldq, 
	float *vl, float *vu, int *il, int *iu, float *abstol, int *
	m, float *w, single_cx *z__, int *ldz, single_cx *work, float *rwork, 
	int *iwork, int *ifail, int *info);
 
int chbgst_(char *vect, char *uplo, int *n, int *ka, 
	int *kb, single_cx *ab, int *ldab, single_cx *bb, int *ldbb, 
	single_cx *x, int *ldx, single_cx *work, float *rwork, int *info);
 
int chbgv_(char *jobz, char *uplo, int *n, int *ka, 
	int *kb, single_cx *ab, int *ldab, single_cx *bb, int *ldbb, 
	float *w, single_cx *z__, int *ldz, single_cx *work, float *rwork, 
	int *info);
 
int chbgvx_(char *jobz, char *range, char *uplo, int *n, 
	int *ka, int *kb, single_cx *ab, int *ldab, single_cx *bb, 
	int *ldbb, single_cx *q, int *ldq, float *vl, float *vu, int *
	il, int *iu, float *abstol, int *m, float *w, single_cx *z__, 
	int *ldz, single_cx *work, float *rwork, int *iwork, int *
	ifail, int *info);
 
int chbtrd_(char *vect, char *uplo, int *n, int *kd, 
	single_cx *ab, int *ldab, float *d__, float *e, single_cx *q, int *
	ldq, single_cx *work, int *info);
 
int checon_(char *uplo, int *n, single_cx *a, int *lda,
	 int *ipiv, float *anorm, float *rcond, single_cx *work, int *
	info);
 
int cheev_(char *jobz, char *uplo, int *n, single_cx *a, 
	int *lda, float *w, single_cx *work, int *lwork, float *rwork, 
	int *info);
 
int cheevd_(char *jobz, char *uplo, int *n, single_cx *a, 
	int *lda, float *w, single_cx *work, int *lwork, float *rwork, 
	int *lrwork, int *iwork, int *liwork, int *info);
 
int cheevr_(char *jobz, char *range, char *uplo, int *n, 
	single_cx *a, int *lda, float *vl, float *vu, int *il, int *
	iu, float *abstol, int *m, float *w, single_cx *z__, int *ldz, 
	int *isuppz, single_cx *work, int *lwork, float *rwork, int *
	lrwork, int *iwork, int *liwork, int *info);
 
int cheevx_(char *jobz, char *range, char *uplo, int *n, 
	single_cx *a, int *lda, float *vl, float *vu, int *il, int *
	iu, float *abstol, int *m, float *w, single_cx *z__, int *ldz, 
	single_cx *work, int *lwork, float *rwork, int *iwork, int *
	ifail, int *info);
 
int chegs2_(int *itype, char *uplo, int *n, single_cx *
	a, int *lda, single_cx *b, int *ldb, int *info);
 
int chegst_(int *itype, char *uplo, int *n, single_cx *
	a, int *lda, single_cx *b, int *ldb, int *info);
 
int chegv_(int *itype, char *jobz, char *uplo, int *
	n, single_cx *a, int *lda, single_cx *b, int *ldb, float *w, 
	single_cx *work, int *lwork, float *rwork, int *info);
 
int chegvd_(int *itype, char *jobz, char *uplo, int *
	n, single_cx *a, int *lda, single_cx *b, int *ldb, float *w, 
	single_cx *work, int *lwork, float *rwork, int *lrwork, int *
	iwork, int *liwork, int *info);
 
int chegvx_(int *itype, char *jobz, char *range, char *
	uplo, int *n, single_cx *a, int *lda, single_cx *b, int *ldb, 
	float *vl, float *vu, int *il, int *iu, float *abstol, int *
	m, float *w, single_cx *z__, int *ldz, single_cx *work, int *lwork,
	 float *rwork, int *iwork, int *ifail, int *info);
 
int cherfs_(char *uplo, int *n, int *nrhs, single_cx *
	a, int *lda, single_cx *af, int *ldaf, int *ipiv, single_cx *
	b, int *ldb, single_cx *x, int *ldx, float *ferr, float *berr, 
	single_cx *work, float *rwork, int *info);
 
int chesv_(char *uplo, int *n, int *nrhs, single_cx *a,
	 int *lda, int *ipiv, single_cx *b, int *ldb, single_cx *work,
	 int *lwork, int *info);
 
int chesvx_(char *fact, char *uplo, int *n, int *
	nrhs, single_cx *a, int *lda, single_cx *af, int *ldaf, int *
	ipiv, single_cx *b, int *ldb, single_cx *x, int *ldx, float *rcond,
	 float *ferr, float *berr, single_cx *work, int *lwork, float *rwork, 
	int *info);
 
int chetf2_(char *uplo, int *n, single_cx *a, int *lda,
	 int *ipiv, int *info);
 
int chetrd_(char *uplo, int *n, single_cx *a, int *lda,
	 float *d__, float *e, single_cx *tau, single_cx *work, int *lwork, 
	int *info);
 
int chetrf_(char *uplo, int *n, single_cx *a, int *lda,
	 int *ipiv, single_cx *work, int *lwork, int *info);
 
int chetri_(char *uplo, int *n, single_cx *a, int *lda,
	 int *ipiv, single_cx *work, int *info);
 
int chetrs_(char *uplo, int *n, int *nrhs, single_cx *
	a, int *lda, int *ipiv, single_cx *b, int *ldb, int *
	info);
 
int chgeqz_(char *job, char *compq, char *compz, int *n, 
	int *ilo, int *ihi, single_cx *a, int *lda, single_cx *b, 
	int *ldb, single_cx *alpha, single_cx *beta, single_cx *q, int *ldq,
	 single_cx *z__, int *ldz, single_cx *work, int *lwork, float *
	rwork, int *info);
 
int chpcon_(char *uplo, int *n, single_cx *ap, int *
	ipiv, float *anorm, float *rcond, single_cx *work, int *info);
 
int chpev_(char *jobz, char *uplo, int *n, single_cx *ap, 
	float *w, single_cx *z__, int *ldz, single_cx *work, float *rwork, 
	int *info);
 
int chpevd_(char *jobz, char *uplo, int *n, single_cx *ap, 
	float *w, single_cx *z__, int *ldz, single_cx *work, int *lwork, 
	float *rwork, int *lrwork, int *iwork, int *liwork, 
	int *info);
 
int chpevx_(char *jobz, char *range, char *uplo, int *n, 
	single_cx *ap, float *vl, float *vu, int *il, int *iu, float *
	abstol, int *m, float *w, single_cx *z__, int *ldz, single_cx *
	work, float *rwork, int *iwork, int *ifail, int *info);
 
int chpgst_(int *itype, char *uplo, int *n, single_cx *
	ap, single_cx *bp, int *info);
 
int chpgv_(int *itype, char *jobz, char *uplo, int *
	n, single_cx *ap, single_cx *bp, float *w, single_cx *z__, int *ldz, 
	single_cx *work, float *rwork, int *info);
 
int chpgvd_(int *itype, char *jobz, char *uplo, int *
	n, single_cx *ap, single_cx *bp, float *w, single_cx *z__, int *ldz, 
	single_cx *work, int *lwork, float *rwork, int *lrwork, int *
	iwork, int *liwork, int *info);
 
int chpgvx_(int *itype, char *jobz, char *range, char *
	uplo, int *n, single_cx *ap, single_cx *bp, float *vl, float *vu, 
	int *il, int *iu, float *abstol, int *m, float *w, single_cx *
	z__, int *ldz, single_cx *work, float *rwork, int *iwork, 
	int *ifail, int *info);
 
int chprfs_(char *uplo, int *n, int *nrhs, single_cx *
	ap, single_cx *afp, int *ipiv, single_cx *b, int *ldb, single_cx *x,
	 int *ldx, float *ferr, float *berr, single_cx *work, float *rwork, 
	int *info);
 
int chpsv_(char *uplo, int *n, int *nrhs, single_cx *
	ap, int *ipiv, single_cx *b, int *ldb, int *info);
 
int chpsvx_(char *fact, char *uplo, int *n, int *
	nrhs, single_cx *ap, single_cx *afp, int *ipiv, single_cx *b, int *
	ldb, single_cx *x, int *ldx, float *rcond, float *ferr, float *berr, 
	single_cx *work, float *rwork, int *info);
 
int chptrd_(char *uplo, int *n, single_cx *ap, float *d__, 
	float *e, single_cx *tau, int *info);
 
int chptrf_(char *uplo, int *n, single_cx *ap, int *
	ipiv, int *info);
 
int chptri_(char *uplo, int *n, single_cx *ap, int *
	ipiv, single_cx *work, int *info);
 
int chptrs_(char *uplo, int *n, int *nrhs, single_cx *
	ap, int *ipiv, single_cx *b, int *ldb, int *info);
 
int chsein_(char *side, char *eigsrc, char *initv, int *
	select, int *n, single_cx *h__, int *ldh, single_cx *w, single_cx *
	vl, int *ldvl, single_cx *vr, int *ldvr, int *mm, int *
	m, single_cx *work, float *rwork, int *ifaill, int *ifailr, 
	int *info);
 
int chseqr_(char *job, char *compz, int *n, int *ilo,
	 int *ihi, single_cx *h__, int *ldh, single_cx *w, single_cx *z__, 
	int *ldz, single_cx *work, int *lwork, int *info);
 
int clabrd_(int *m, int *n, int *nb, single_cx *a, 
	int *lda, float *d__, float *e, single_cx *tauq, single_cx *taup, 
	single_cx *x, int *ldx, single_cx *y, int *ldy);
 
int clacgv_(int *n, single_cx *x, int *incx);
 
int clacon_(int *n, single_cx *v, single_cx *x, float *est, 
	int *kase);
 
int clacp2_(char *uplo, int *m, int *n, float *a, 
	int *lda, single_cx *b, int *ldb);
 
int clacpy_(char *uplo, int *m, int *n, single_cx *a, 
	int *lda, single_cx *b, int *ldb);
 
int clacrm_(int *m, int *n, single_cx *a, int *lda,
	 float *b, int *ldb, single_cx *c__, int *ldc, float *rwork);
 
int clacrt_(int *n, single_cx *cx, int *incx, single_cx *
	cy, int *incy, single_cx *c__, single_cx *s);
 
int claed0_(int *qsiz, int *n, float *d__, float *e, 
	single_cx *q, int *ldq, single_cx *qstore, int *ldqs, float *rwork,
	 int *iwork, int *info);
 
int claed7_(int *n, int *cutpnt, int *qsiz, 
	int *tlvls, int *curlvl, int *curpbm, float *d__, single_cx *
	q, int *ldq, float *rho, int *indxq, float *qstore, int *
	qptr, int *prmptr, int *perm, int *givptr, int *
	givcol, float *givnum, single_cx *work, float *rwork, int *iwork, 
	int *info);
 
int claed8_(int *k, int *n, int *qsiz, single_cx *
	q, int *ldq, float *d__, float *rho, int *cutpnt, float *z__, 
	float *dlamda, single_cx *q2, int *ldq2, float *w, int *indxp, 
	int *indx, int *indxq, int *perm, int *givptr, 
	int *givcol, float *givnum, int *info);
 
int claein_(int *rightv, int *noinit, int *n, 
	single_cx *h__, int *ldh, single_cx *w, single_cx *v, single_cx *b, 
	int *ldb, float *rwork, float *eps3, float *smlnum, int *info);
 
int claesy_(single_cx *a, single_cx *b, single_cx *c__, single_cx *
	rt1, single_cx *rt2, single_cx *evscal, single_cx *cs1, single_cx *sn1);
 
int claev2_(single_cx *a, single_cx *b, single_cx *c__, float *rt1, 
	float *rt2, float *cs1, single_cx *sn1);
 
int clags2_(int *upper, float *a1, single_cx *a2, float *a3, 
	float *b1, single_cx *b2, float *b3, float *csu, single_cx *snu, float *csv, 
	single_cx *snv, float *csq, single_cx *snq);
 
int clagtm_(char *trans, int *n, int *nrhs, float *
	alpha, single_cx *dl, single_cx *d__, single_cx *du, single_cx *x, int *
	ldx, float *beta, single_cx *b, int *ldb);
 
int clahef_(char *uplo, int *n, int *nb, int *kb,
	 single_cx *a, int *lda, int *ipiv, single_cx *w, int *ldw, 
	int *info);
 
int clahqr_(int *wantt, int *wantz, int *n, 
	int *ilo, int *ihi, single_cx *h__, int *ldh, single_cx *w, 
	int *iloz, int *ihiz, single_cx *z__, int *ldz, int *
	info);
 
int clahrd_(int *n, int *k, int *nb, single_cx *a, 
	int *lda, single_cx *tau, single_cx *t, int *ldt, single_cx *y, 
	int *ldy);
 
int claic1_(int *job, int *j, single_cx *x, float *sest,
	 single_cx *w, single_cx *gamma, float *sestpr, single_cx *s, single_cx *c__);
 
int clals0_(int *icompq, int *nl, int *nr, 
	int *sqre, int *nrhs, single_cx *b, int *ldb, single_cx *bx, 
	int *ldbx, int *perm, int *givptr, int *givcol, 
	int *ldgcol, float *givnum, int *ldgnum, float *poles, float *
	difl, float *difr, float *z__, int *k, float *c__, float *s, float *
	rwork, int *info);
 
int clalsa_(int *icompq, int *smlsiz, int *n, 
	int *nrhs, single_cx *b, int *ldb, single_cx *bx, int *ldbx, 
	float *u, int *ldu, float *vt, int *k, float *difl, float *difr, 
	float *z__, float *poles, int *givptr, int *givcol, int *
	ldgcol, int *perm, float *givnum, float *c__, float *s, float *rwork, 
	int *iwork, int *info);
 
int clapll_(int *n, single_cx *x, int *incx, single_cx *
	y, int *incy, float *ssmin);
 
int clapmt_(int *forwrd, int *m, int *n, single_cx 
	*x, int *ldx, int *k);
 
int claqgb_(int *m, int *n, int *kl, int *ku,
	 single_cx *ab, int *ldab, float *r__, float *c__, float *rowcnd, float 
	*colcnd, float *amax, char *equed);
 
int claqge_(int *m, int *n, single_cx *a, int *lda,
	 float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, char *
	equed);
 
int claqhb_(char *uplo, int *n, int *kd, single_cx *ab,
	 int *ldab, float *s, float *scond, float *amax, char *equed);
 
int claqhe_(char *uplo, int *n, single_cx *a, int *lda,
	 float *s, float *scond, float *amax, char *equed);
 
int claqhp_(char *uplo, int *n, single_cx *ap, float *s, 
	float *scond, float *amax, char *equed);
 
int claqp2_(int *m, int *n, int *offset, single_cx 
	*a, int *lda, int *jpvt, single_cx *tau, float *vn1, float *vn2, 
	single_cx *work);
 
int claqps_(int *m, int *n, int *offset, int 
	*nb, int *kb, single_cx *a, int *lda, int *jpvt, single_cx *
	tau, float *vn1, float *vn2, single_cx *auxv, single_cx *f, int *ldf);
 
int claqsb_(char *uplo, int *n, int *kd, single_cx *ab,
	 int *ldab, float *s, float *scond, float *amax, char *equed);
 
int claqsp_(char *uplo, int *n, single_cx *ap, float *s, 
	float *scond, float *amax, char *equed);
 
int claqsy_(char *uplo, int *n, single_cx *a, int *lda,
	 float *s, float *scond, float *amax, char *equed);
 
int clar1v_(int *n, int *b1, int *bn, float *
	sigma, float *d__, float *l, float *ld, float *lld, float *gersch, single_cx 
	*z__, float *ztz, float *mingma, int *r__, int *isuppz, float *
	work);
 
int clar2v_(int *n, single_cx *x, single_cx *y, single_cx *z__,
	 int *incx, float *c__, single_cx *s, int *incc);
 
int clarcm_(int *m, int *n, float *a, int *lda, 
	single_cx *b, int *ldb, single_cx *c__, int *ldc, float *rwork);
 
int clarf_(char *side, int *m, int *n, single_cx *v, 
	int *incv, single_cx *tau, single_cx *c__, int *ldc, single_cx *
	work);
 
int clarfb_(char *side, char *trans, char *direct, char *
	storev, int *m, int *n, int *k, single_cx *v, int *ldv, 
	single_cx *t, int *ldt, single_cx *c__, int *ldc, single_cx *work, 
	int *ldwork);
 
int clarfg_(int *n, single_cx *alpha, single_cx *x, int *
	incx, single_cx *tau);
 
int clarft_(char *direct, char *storev, int *n, int *
	k, single_cx *v, int *ldv, single_cx *tau, single_cx *t, int *ldt);
 
int clarfx_(char *side, int *m, int *n, single_cx *v, 
	single_cx *tau, single_cx *c__, int *ldc, single_cx *work);
 
int clargv_(int *n, single_cx *x, int *incx, single_cx *
	y, int *incy, float *c__, int *incc);
 
int clarnv_(int *idist, int *iseed, int *n, 
	single_cx *x);
 
int clarrv_(int *n, float *d__, float *l, int *isplit, 
	int *m, float *w, int *iblock, float *gersch, float *tol, 
	single_cx *z__, int *ldz, int *isuppz, float *work, int *
	iwork, int *info);
 
int clartg_(single_cx *f, single_cx *g, float *cs, single_cx *sn, 
	single_cx *r__);
 
int clartv_(int *n, single_cx *x, int *incx, single_cx *
	y, int *incy, float *c__, single_cx *s, int *incc);
 
int clarz_(char *side, int *m, int *n, int *l, 
	single_cx *v, int *incv, single_cx *tau, single_cx *c__, int *ldc, 
	single_cx *work);
 
int clarzb_(char *side, char *trans, char *direct, char *
	storev, int *m, int *n, int *k, int *l, single_cx *v, 
	int *ldv, single_cx *t, int *ldt, single_cx *c__, int *ldc, 
	single_cx *work, int *ldwork);
 
int clarzt_(char *direct, char *storev, int *n, int *
	k, single_cx *v, int *ldv, single_cx *tau, single_cx *t, int *ldt);
 
int clascl_(char *type__, int *kl, int *ku, float *
	cfrom, float *cto, int *m, int *n, single_cx *a, int *lda, 
	int *info);
 
int claset_(char *uplo, int *m, int *n, single_cx *
	alpha, single_cx *beta, single_cx *a, int *lda);
 
int clasr_(char *side, char *pivot, char *direct, int *m,
	 int *n, float *c__, float *s, single_cx *a, int *lda);
 
int classq_(int *n, single_cx *x, int *incx, float *
	scale, float *sumsq);
 
int claswp_(int *n, single_cx *a, int *lda, int *
	k1, int *k2, int *ipiv, int *incx);
 
int clasyf_(char *uplo, int *n, int *nb, int *kb,
	 single_cx *a, int *lda, int *ipiv, single_cx *w, int *ldw, 
	int *info);
 
int clatbs_(char *uplo, char *trans, char *diag, char *
	normin, int *n, int *kd, single_cx *ab, int *ldab, single_cx *
	x, float *scale, float *cnorm, int *info);
 
int clatdf_(int *ijob, int *n, single_cx *z__, int 
	*ldz, single_cx *rhs, float *rdsum, float *rdscal, int *ipiv, int 
	*jpiv);
 
int clatps_(char *uplo, char *trans, char *diag, char *
	normin, int *n, single_cx *ap, single_cx *x, float *scale, float *cnorm,
	 int *info);
 
int clatrd_(char *uplo, int *n, int *nb, single_cx *a, 
	int *lda, float *e, single_cx *tau, single_cx *w, int *ldw);
 
int clatrs_(char *uplo, char *trans, char *diag, char *
	normin, int *n, single_cx *a, int *lda, single_cx *x, float *scale,
	 float *cnorm, int *info);
 
int clatrz_(int *m, int *n, int *l, single_cx *a, 
	int *lda, single_cx *tau, single_cx *work);
 
int clatzm_(char *side, int *m, int *n, single_cx *v, 
	int *incv, single_cx *tau, single_cx *c1, single_cx *c2, int *ldc, 
	single_cx *work);
 
int clauu2_(char *uplo, int *n, single_cx *a, int *lda,
	 int *info);
 
int clauum_(char *uplo, int *n, single_cx *a, int *lda,
	 int *info);
 
int cpbcon_(char *uplo, int *n, int *kd, single_cx *ab,
	 int *ldab, float *anorm, float *rcond, single_cx *work, float *rwork, 
	int *info);
 
int cpbequ_(char *uplo, int *n, int *kd, single_cx *ab,
	 int *ldab, float *s, float *scond, float *amax, int *info);
 
int cpbrfs_(char *uplo, int *n, int *kd, int *
	nrhs, single_cx *ab, int *ldab, single_cx *afb, int *ldafb, 
	single_cx *b, int *ldb, single_cx *x, int *ldx, float *ferr, float *
	berr, single_cx *work, float *rwork, int *info);
 
int cpbstf_(char *uplo, int *n, int *kd, single_cx *ab,
	 int *ldab, int *info);
 
int cpbsv_(char *uplo, int *n, int *kd, int *
	nrhs, single_cx *ab, int *ldab, single_cx *b, int *ldb, int *
	info);
 
int cpbsvx_(char *fact, char *uplo, int *n, int *kd, 
	int *nrhs, single_cx *ab, int *ldab, single_cx *afb, int *
	ldafb, char *equed, float *s, single_cx *b, int *ldb, single_cx *x, 
	int *ldx, float *rcond, float *ferr, float *berr, single_cx *work, 
	float *rwork, int *info);
 
int cpbtf2_(char *uplo, int *n, int *kd, single_cx *ab,
	 int *ldab, int *info);
 
int cpbtrf_(char *uplo, int *n, int *kd, single_cx *ab,
	 int *ldab, int *info);
 
int cpbtrs_(char *uplo, int *n, int *kd, int *
	nrhs, single_cx *ab, int *ldab, single_cx *b, int *ldb, int *
	info);
 
int cpocon_(char *uplo, int *n, single_cx *a, int *lda,
	 float *anorm, float *rcond, single_cx *work, float *rwork, int *info);
 
int cpoequ_(int *n, single_cx *a, int *lda, float *s, 
	float *scond, float *amax, int *info);
 
int cporfs_(char *uplo, int *n, int *nrhs, single_cx *
	a, int *lda, single_cx *af, int *ldaf, single_cx *b, int *ldb,
	 single_cx *x, int *ldx, float *ferr, float *berr, single_cx *work, 
	float *rwork, int *info);
 
int cposv_(char *uplo, int *n, int *nrhs, single_cx *a,
	 int *lda, single_cx *b, int *ldb, int *info);
 
int cposvx_(char *fact, char *uplo, int *n, int *
	nrhs, single_cx *a, int *lda, single_cx *af, int *ldaf, char *
	equed, float *s, single_cx *b, int *ldb, single_cx *x, int *ldx, 
	float *rcond, float *ferr, float *berr, single_cx *work, float *rwork, 
	int *info);
 
int cpotf2_(char *uplo, int *n, single_cx *a, int *lda,
	 int *info);
 
int cpotrf_(char *uplo, int *n, single_cx *a, int *lda,
	 int *info);
 
int cpotri_(char *uplo, int *n, single_cx *a, int *lda,
	 int *info);
 
int cpotrs_(char *uplo, int *n, int *nrhs, single_cx *
	a, int *lda, single_cx *b, int *ldb, int *info);
 
int cppcon_(char *uplo, int *n, single_cx *ap, float *anorm,
	 float *rcond, single_cx *work, float *rwork, int *info);
 
int cppequ_(char *uplo, int *n, single_cx *ap, float *s, 
	float *scond, float *amax, int *info);
 
int cpprfs_(char *uplo, int *n, int *nrhs, single_cx *
	ap, single_cx *afp, single_cx *b, int *ldb, single_cx *x, int *ldx, 
	float *ferr, float *berr, single_cx *work, float *rwork, int *info);
 
int cppsv_(char *uplo, int *n, int *nrhs, single_cx *
	ap, single_cx *b, int *ldb, int *info);
 
int cppsvx_(char *fact, char *uplo, int *n, int *
	nrhs, single_cx *ap, single_cx *afp, char *equed, float *s, single_cx *b, 
	int *ldb, single_cx *x, int *ldx, float *rcond, float *ferr, float 
	*berr, single_cx *work, float *rwork, int *info);
 
int cpptrf_(char *uplo, int *n, single_cx *ap, int *
	info);
 
int cpptri_(char *uplo, int *n, single_cx *ap, int *
	info);
 
int cpptrs_(char *uplo, int *n, int *nrhs, single_cx *
	ap, single_cx *b, int *ldb, int *info);
 
int cptcon_(int *n, float *d__, single_cx *e, float *anorm, 
	float *rcond, float *rwork, int *info);
 
int cptrfs_(char *uplo, int *n, int *nrhs, float *d__,
	 single_cx *e, float *df, single_cx *ef, single_cx *b, int *ldb, single_cx 
	*x, int *ldx, float *ferr, float *berr, single_cx *work, float *rwork, 
	int *info);
 
int cptsv_(int *n, int *nrhs, float *d__, single_cx *e, 
	single_cx *b, int *ldb, int *info);
 
int cptsvx_(char *fact, int *n, int *nrhs, float *d__,
	 single_cx *e, float *df, single_cx *ef, single_cx *b, int *ldb, single_cx 
	*x, int *ldx, float *rcond, float *ferr, float *berr, single_cx *work, 
	float *rwork, int *info);
 
int cpttrf_(int *n, float *d__, single_cx *e, int *info);
 
int cpttrs_(char *uplo, int *n, int *nrhs, float *d__,
	 single_cx *e, single_cx *b, int *ldb, int *info);
 
int cptts2_(int *iuplo, int *n, int *nrhs, float *
	d__, single_cx *e, single_cx *b, int *ldb);
 
int crot_(int *n, single_cx *cx, int *incx, single_cx *
	cy, int *incy, float *c__, single_cx *s);
 
int cspcon_(char *uplo, int *n, single_cx *ap, int *
	ipiv, float *anorm, float *rcond, single_cx *work, int *info);
 
int cspmv_(char *uplo, int *n, single_cx *alpha, single_cx *
	ap, single_cx *x, int *incx, single_cx *beta, single_cx *y, int *
	incy);
 
int cspr_(char *uplo, int *n, single_cx *alpha, single_cx *x,
	 int *incx, single_cx *ap);
 
int csprfs_(char *uplo, int *n, int *nrhs, single_cx *
	ap, single_cx *afp, int *ipiv, single_cx *b, int *ldb, single_cx *x,
	 int *ldx, float *ferr, float *berr, single_cx *work, float *rwork, 
	int *info);
 
int cspsv_(char *uplo, int *n, int *nrhs, single_cx *
	ap, int *ipiv, single_cx *b, int *ldb, int *info);
 
int cspsvx_(char *fact, char *uplo, int *n, int *
	nrhs, single_cx *ap, single_cx *afp, int *ipiv, single_cx *b, int *
	ldb, single_cx *x, int *ldx, float *rcond, float *ferr, float *berr, 
	single_cx *work, float *rwork, int *info);
 
int csptrf_(char *uplo, int *n, single_cx *ap, int *
	ipiv, int *info);
 
int csptri_(char *uplo, int *n, single_cx *ap, int *
	ipiv, single_cx *work, int *info);
 
int csptrs_(char *uplo, int *n, int *nrhs, single_cx *
	ap, int *ipiv, single_cx *b, int *ldb, int *info);
 
int csrot_(int *n, single_cx *cx, int *incx, single_cx *
	cy, int *incy, float *c__, float *s);
 
int csrscl_(int *n, float *sa, single_cx *sx, int *incx);
 
int cstedc_(char *compz, int *n, float *d__, float *e, 
	single_cx *z__, int *ldz, single_cx *work, int *lwork, float *
	rwork, int *lrwork, int *iwork, int *liwork, int *
	info);
 
int cstein_(int *n, float *d__, float *e, int *m, float 
	*w, int *iblock, int *isplit, single_cx *z__, int *ldz, 
	float *work, int *iwork, int *ifail, int *info);
 
int csteqr_(char *compz, int *n, float *d__, float *e, 
	single_cx *z__, int *ldz, float *work, int *info);
 
int csycon_(char *uplo, int *n, single_cx *a, int *lda,
	 int *ipiv, float *anorm, float *rcond, single_cx *work, int *
	info);
 
int csymv_(char *uplo, int *n, single_cx *alpha, single_cx *
	a, int *lda, single_cx *x, int *incx, single_cx *beta, single_cx *y,
	 int *incy);
 
int csyr_(char *uplo, int *n, single_cx *alpha, single_cx *x,
	 int *incx, single_cx *a, int *lda);
 
int csyrfs_(char *uplo, int *n, int *nrhs, single_cx *
	a, int *lda, single_cx *af, int *ldaf, int *ipiv, single_cx *
	b, int *ldb, single_cx *x, int *ldx, float *ferr, float *berr, 
	single_cx *work, float *rwork, int *info);
 
int csysv_(char *uplo, int *n, int *nrhs, single_cx *a,
	 int *lda, int *ipiv, single_cx *b, int *ldb, single_cx *work,
	 int *lwork, int *info);
 
int csysvx_(char *fact, char *uplo, int *n, int *
	nrhs, single_cx *a, int *lda, single_cx *af, int *ldaf, int *
	ipiv, single_cx *b, int *ldb, single_cx *x, int *ldx, float *rcond,
	 float *ferr, float *berr, single_cx *work, int *lwork, float *rwork, 
	int *info);
 
int csytf2_(char *uplo, int *n, single_cx *a, int *lda,
	 int *ipiv, int *info);
 
int csytrf_(char *uplo, int *n, single_cx *a, int *lda,
	 int *ipiv, single_cx *work, int *lwork, int *info);
 
int csytri_(char *uplo, int *n, single_cx *a, int *lda,
	 int *ipiv, single_cx *work, int *info);
 
int csytrs_(char *uplo, int *n, int *nrhs, single_cx *
	a, int *lda, int *ipiv, single_cx *b, int *ldb, int *
	info);
 
int ctbcon_(char *norm, char *uplo, char *diag, int *n, 
	int *kd, single_cx *ab, int *ldab, float *rcond, single_cx *work, 
	float *rwork, int *info);
 
int ctbrfs_(char *uplo, char *trans, char *diag, int *n, 
	int *kd, int *nrhs, single_cx *ab, int *ldab, single_cx *b, 
	int *ldb, single_cx *x, int *ldx, float *ferr, float *berr, 
	single_cx *work, float *rwork, int *info);
 
int ctbtrs_(char *uplo, char *trans, char *diag, int *n, 
	int *kd, int *nrhs, single_cx *ab, int *ldab, single_cx *b, 
	int *ldb, int *info);
 
int ctgevc_(char *side, char *howmny, int *select, 
	int *n, single_cx *a, int *lda, single_cx *b, int *ldb, 
	single_cx *vl, int *ldvl, single_cx *vr, int *ldvr, int *mm, 
	int *m, single_cx *work, float *rwork, int *info);
 
int ctgex2_(int *wantq, int *wantz, int *n, 
	single_cx *a, int *lda, single_cx *b, int *ldb, single_cx *q, 
	int *ldq, single_cx *z__, int *ldz, int *j1, int *info);
 
int ctgexc_(int *wantq, int *wantz, int *n, 
	single_cx *a, int *lda, single_cx *b, int *ldb, single_cx *q, 
	int *ldq, single_cx *z__, int *ldz, int *ifst, int *
	ilst, int *info);
 
int ctgsen_(int *ijob, int *wantq, int *wantz, 
	int *select, int *n, single_cx *a, int *lda, single_cx *b, 
	int *ldb, single_cx *alpha, single_cx *beta, single_cx *q, int *ldq,
	 single_cx *z__, int *ldz, int *m, float *pl, float *pr, float *
	dif, single_cx *work, int *lwork, int *iwork, int *liwork, 
	int *info);
 
int ctgsja_(char *jobu, char *jobv, char *jobq, int *m, 
	int *p, int *n, int *k, int *l, single_cx *a, int *
	lda, single_cx *b, int *ldb, float *tola, float *tolb, float *alpha, 
	float *beta, single_cx *u, int *ldu, single_cx *v, int *ldv, 
	single_cx *q, int *ldq, single_cx *work, int *ncycle, int *
	info);
 
int ctgsna_(char *job, char *howmny, int *select, 
	int *n, single_cx *a, int *lda, single_cx *b, int *ldb, 
	single_cx *vl, int *ldvl, single_cx *vr, int *ldvr, float *s, float 
	*dif, int *mm, int *m, single_cx *work, int *lwork, int 
	*iwork, int *info);
 
int ctgsy2_(char *trans, int *ijob, int *m, int *
	n, single_cx *a, int *lda, single_cx *b, int *ldb, single_cx *c__, 
	int *ldc, single_cx *d__, int *ldd, single_cx *e, int *lde, 
	single_cx *f, int *ldf, float *scale, float *rdsum, float *rdscal, 
	int *info);
 
int ctgsyl_(char *trans, int *ijob, int *m, int *
	n, single_cx *a, int *lda, single_cx *b, int *ldb, single_cx *c__, 
	int *ldc, single_cx *d__, int *ldd, single_cx *e, int *lde, 
	single_cx *f, int *ldf, float *scale, float *dif, single_cx *work, 
	int *lwork, int *iwork, int *info);
 
int ctpcon_(char *norm, char *uplo, char *diag, int *n, 
	single_cx *ap, float *rcond, single_cx *work, float *rwork, int *info);
 
int ctprfs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, single_cx *ap, single_cx *b, int *ldb, single_cx *x, 
	int *ldx, float *ferr, float *berr, single_cx *work, float *rwork, 
	int *info);
 
int ctptri_(char *uplo, char *diag, int *n, single_cx *ap, 
	int *info);
 
int ctptrs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, single_cx *ap, single_cx *b, int *ldb, int *info);
 
int ctrcon_(char *norm, char *uplo, char *diag, int *n, 
	single_cx *a, int *lda, float *rcond, single_cx *work, float *rwork, 
	int *info);
 
int ctrevc_(char *side, char *howmny, int *select, 
	int *n, single_cx *t, int *ldt, single_cx *vl, int *ldvl, 
	single_cx *vr, int *ldvr, int *mm, int *m, single_cx *work, 
	float *rwork, int *info);
 
int ctrexc_(char *compq, int *n, single_cx *t, int *
	ldt, single_cx *q, int *ldq, int *ifst, int *ilst, int *
	info);
 
int ctrrfs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, single_cx *a, int *lda, single_cx *b, int *ldb, 
	single_cx *x, int *ldx, float *ferr, float *berr, single_cx *work, float 
	*rwork, int *info);
 
int ctrsen_(char *job, char *compq, int *select, int 
	*n, single_cx *t, int *ldt, single_cx *q, int *ldq, single_cx *w, 
	int *m, float *s, float *sep, single_cx *work, int *lwork, 
	int *info);
 
int ctrsna_(char *job, char *howmny, int *select, 
	int *n, single_cx *t, int *ldt, single_cx *vl, int *ldvl, 
	single_cx *vr, int *ldvr, float *s, float *sep, int *mm, int *
	m, single_cx *work, int *ldwork, float *rwork, int *info);
 
int ctrsyl_(char *trana, char *tranb, int *isgn, int 
	*m, int *n, single_cx *a, int *lda, single_cx *b, int *ldb, 
	single_cx *c__, int *ldc, float *scale, int *info);
 
int ctrti2_(char *uplo, char *diag, int *n, single_cx *a, 
	int *lda, int *info);
 
int ctrtri_(char *uplo, char *diag, int *n, single_cx *a, 
	int *lda, int *info);
 
int ctrtrs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, single_cx *a, int *lda, single_cx *b, int *ldb, 
	int *info);
 
int ctzrqf_(int *m, int *n, single_cx *a, int *lda,
	 single_cx *tau, int *info);
 
int ctzrzf_(int *m, int *n, single_cx *a, int *lda,
	 single_cx *tau, single_cx *work, int *lwork, int *info);
 
int cung2l_(int *m, int *n, int *k, single_cx *a, 
	int *lda, single_cx *tau, single_cx *work, int *info);
 
int cung2r_(int *m, int *n, int *k, single_cx *a, 
	int *lda, single_cx *tau, single_cx *work, int *info);
 
int cungbr_(char *vect, int *m, int *n, int *k, 
	single_cx *a, int *lda, single_cx *tau, single_cx *work, int *lwork,
	 int *info);
 
int cunghr_(int *n, int *ilo, int *ihi, single_cx *
	a, int *lda, single_cx *tau, single_cx *work, int *lwork, int 
	*info);
 
int cungl2_(int *m, int *n, int *k, single_cx *a, 
	int *lda, single_cx *tau, single_cx *work, int *info);
 
int cunglq_(int *m, int *n, int *k, single_cx *a, 
	int *lda, single_cx *tau, single_cx *work, int *lwork, int *
	info);
 
int cungql_(int *m, int *n, int *k, single_cx *a, 
	int *lda, single_cx *tau, single_cx *work, int *lwork, int *
	info);
 
int cungqr_(int *m, int *n, int *k, single_cx *a, 
	int *lda, single_cx *tau, single_cx *work, int *lwork, int *
	info);
 
int cungr2_(int *m, int *n, int *k, single_cx *a, 
	int *lda, single_cx *tau, single_cx *work, int *info);
 
int cungrq_(int *m, int *n, int *k, single_cx *a, 
	int *lda, single_cx *tau, single_cx *work, int *lwork, int *
	info);
 
int cungtr_(char *uplo, int *n, single_cx *a, int *lda,
	 single_cx *tau, single_cx *work, int *lwork, int *info);
 
int cunm2l_(char *side, char *trans, int *m, int *n, 
	int *k, single_cx *a, int *lda, single_cx *tau, single_cx *c__, 
	int *ldc, single_cx *work, int *info);
 
int cunm2r_(char *side, char *trans, int *m, int *n, 
	int *k, single_cx *a, int *lda, single_cx *tau, single_cx *c__, 
	int *ldc, single_cx *work, int *info);
 
int cunmbr_(char *vect, char *side, char *trans, int *m, 
	int *n, int *k, single_cx *a, int *lda, single_cx *tau, 
	single_cx *c__, int *ldc, single_cx *work, int *lwork, int *
	info);
 
int cunmhr_(char *side, char *trans, int *m, int *n, 
	int *ilo, int *ihi, single_cx *a, int *lda, single_cx *tau, 
	single_cx *c__, int *ldc, single_cx *work, int *lwork, int *
	info);
 
int cunml2_(char *side, char *trans, int *m, int *n, 
	int *k, single_cx *a, int *lda, single_cx *tau, single_cx *c__, 
	int *ldc, single_cx *work, int *info);
 
int cunmlq_(char *side, char *trans, int *m, int *n, 
	int *k, single_cx *a, int *lda, single_cx *tau, single_cx *c__, 
	int *ldc, single_cx *work, int *lwork, int *info);
 
int cunmql_(char *side, char *trans, int *m, int *n, 
	int *k, single_cx *a, int *lda, single_cx *tau, single_cx *c__, 
	int *ldc, single_cx *work, int *lwork, int *info);
 
int cunmqr_(char *side, char *trans, int *m, int *n, 
	int *k, single_cx *a, int *lda, single_cx *tau, single_cx *c__, 
	int *ldc, single_cx *work, int *lwork, int *info);
 
int cunmr2_(char *side, char *trans, int *m, int *n, 
	int *k, single_cx *a, int *lda, single_cx *tau, single_cx *c__, 
	int *ldc, single_cx *work, int *info);
 
int cunmr3_(char *side, char *trans, int *m, int *n, 
	int *k, int *l, single_cx *a, int *lda, single_cx *tau, 
	single_cx *c__, int *ldc, single_cx *work, int *info);
 
int cunmrq_(char *side, char *trans, int *m, int *n, 
	int *k, single_cx *a, int *lda, single_cx *tau, single_cx *c__, 
	int *ldc, single_cx *work, int *lwork, int *info);
 
int cunmrz_(char *side, char *trans, int *m, int *n, 
	int *k, int *l, single_cx *a, int *lda, single_cx *tau, 
	single_cx *c__, int *ldc, single_cx *work, int *lwork, int *
	info);
 
int cunmtr_(char *side, char *uplo, char *trans, int *m, 
	int *n, single_cx *a, int *lda, single_cx *tau, single_cx *c__, 
	int *ldc, single_cx *work, int *lwork, int *info);
 
int cupgtr_(char *uplo, int *n, single_cx *ap, single_cx *
	tau, single_cx *q, int *ldq, single_cx *work, int *info);
 
int cupmtr_(char *side, char *uplo, char *trans, int *m, 
	int *n, single_cx *ap, single_cx *tau, single_cx *c__, int *ldc, 
	single_cx *work, int *info);
 
int dbdsdc_(char *uplo, char *compq, int *n, double *
	d__, double *e, double *u, int *ldu, double *vt, 
	int *ldvt, double *q, int *iq, double *work, int *
	iwork, int *info);
 
int dbdsqr_(char *uplo, int *n, int *ncvt, int *
	nru, int *ncc, double *d__, double *e, double *vt, 
	int *ldvt, double *u, int *ldu, double *c__, int *
	ldc, double *work, int *info);
 
int ddisna_(char *job, int *m, int *n, double *
	d__, double *sep, int *info);
 
int dgbbrd_(char *vect, int *m, int *n, int *ncc,
	 int *kl, int *ku, double *ab, int *ldab, double *
	d__, double *e, double *q, int *ldq, double *pt, 
	int *ldpt, double *c__, int *ldc, double *work, 
	int *info);
 
int dgbcon_(char *norm, int *n, int *kl, int *ku,
	 double *ab, int *ldab, int *ipiv, double *anorm, 
	double *rcond, double *work, int *iwork, int *info);
 
int dgbequ_(int *m, int *n, int *kl, int *ku,
	 double *ab, int *ldab, double *r__, double *c__, 
	double *rowcnd, double *colcnd, double *amax, int *
	info);
 
int dgbrfs_(char *trans, int *n, int *kl, int *
	ku, int *nrhs, double *ab, int *ldab, double *afb, 
	int *ldafb, int *ipiv, double *b, int *ldb, 
	double *x, int *ldx, double *ferr, double *berr, 
	double *work, int *iwork, int *info);
 
int dgbsv_(int *n, int *kl, int *ku, int *
	nrhs, double *ab, int *ldab, int *ipiv, double *b, 
	int *ldb, int *info);
 
int dgbsvx_(char *fact, char *trans, int *n, int *kl,
	 int *ku, int *nrhs, double *ab, int *ldab, 
	double *afb, int *ldafb, int *ipiv, char *equed, 
	double *r__, double *c__, double *b, int *ldb, 
	double *x, int *ldx, double *rcond, double *ferr, 
	double *berr, double *work, int *iwork, int *info);
 
int dgbtf2_(int *m, int *n, int *kl, int *ku,
	 double *ab, int *ldab, int *ipiv, int *info);
 
int dgbtrf_(int *m, int *n, int *kl, int *ku,
	 double *ab, int *ldab, int *ipiv, int *info);
 
int dgbtrs_(char *trans, int *n, int *kl, int *
	ku, int *nrhs, double *ab, int *ldab, int *ipiv, 
	double *b, int *ldb, int *info);
 
int dgebak_(char *job, char *side, int *n, int *ilo, 
	int *ihi, double *scale, int *m, double *v, int *
	ldv, int *info);
 
int dgebal_(char *job, int *n, double *a, int *
	lda, int *ilo, int *ihi, double *scale, int *info);
 
int dgebd2_(int *m, int *n, double *a, int *
	lda, double *d__, double *e, double *tauq, double *
	taup, double *work, int *info);
 
int dgebrd_(int *m, int *n, double *a, int *
	lda, double *d__, double *e, double *tauq, double *
	taup, double *work, int *lwork, int *info);
 
int dgecon_(char *norm, int *n, double *a, int *
	lda, double *anorm, double *rcond, double *work, int *
	iwork, int *info);
 
int dgeequ_(int *m, int *n, double *a, int *
	lda, double *r__, double *c__, double *rowcnd, double 
	*colcnd, double *amax, int *info);
 
int dgees_(char *jobvs, char *sort, L_fp select, int *n, 
	double *a, int *lda, int *sdim, double *wr, 
	double *wi, double *vs, int *ldvs, double *work, 
	int *lwork, int *bwork, int *info);
 
int dgeesx_(char *jobvs, char *sort, L_fp select, char *
	sense, int *n, double *a, int *lda, int *sdim, 
	double *wr, double *wi, double *vs, int *ldvs, 
	double *rconde, double *rcondv, double *work, int *
	lwork, int *iwork, int *liwork, int *bwork, int *info);
 
int dgeev_(char *jobvl, char *jobvr, int *n, double *
	a, int *lda, double *wr, double *wi, double *vl, 
	int *ldvl, double *vr, int *ldvr, double *work, 
	int *lwork, int *info);
 
int dgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, int *n, double *a, int *lda, double *wr, 
	double *wi, double *vl, int *ldvl, double *vr, 
	int *ldvr, int *ilo, int *ihi, double *scale, 
	double *abnrm, double *rconde, double *rcondv, double 
	*work, int *lwork, int *iwork, int *info);
 
int dgegs_(char *jobvsl, char *jobvsr, int *n, 
	double *a, int *lda, double *b, int *ldb, double *
	alphar, double *alphai, double *beta, double *vsl, 
	int *ldvsl, double *vsr, int *ldvsr, double *work, 
	int *lwork, int *info);
 
int dgegv_(char *jobvl, char *jobvr, int *n, double *
	a, int *lda, double *b, int *ldb, double *alphar, 
	double *alphai, double *beta, double *vl, int *ldvl, 
	double *vr, int *ldvr, double *work, int *lwork, 
	int *info);
 
int dgehd2_(int *n, int *ilo, int *ihi, 
	double *a, int *lda, double *tau, double *work, 
	int *info);
 
int dgehrd_(int *n, int *ilo, int *ihi, 
	double *a, int *lda, double *tau, double *work, 
	int *lwork, int *info);
 
int dgelq2_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *info);
 
int dgelqf_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *lwork, int *info);
 
int dgels_(char *trans, int *m, int *n, int *
	nrhs, double *a, int *lda, double *b, int *ldb, 
	double *work, int *lwork, int *info);
 
int dgelsd_(int *m, int *n, int *nrhs, 
	double *a, int *lda, double *b, int *ldb, double *
	s, double *rcond, int *rank, double *work, int *lwork,
	 int *iwork, int *info);
 
int dgelss_(int *m, int *n, int *nrhs, 
	double *a, int *lda, double *b, int *ldb, double *
	s, double *rcond, int *rank, double *work, int *lwork,
	 int *info);
 
int dgelsx_(int *m, int *n, int *nrhs, 
	double *a, int *lda, double *b, int *ldb, int *
	jpvt, double *rcond, int *rank, double *work, int *
	info);
 
int dgelsy_(int *m, int *n, int *nrhs, 
	double *a, int *lda, double *b, int *ldb, int *
	jpvt, double *rcond, int *rank, double *work, int *
	lwork, int *info);
 
int dgeql2_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *info);
 
int dgeqlf_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *lwork, int *info);
 
int dgeqp3_(int *m, int *n, double *a, int *
	lda, int *jpvt, double *tau, double *work, int *lwork,
	 int *info);
 
int dgeqpf_(int *m, int *n, double *a, int *
	lda, int *jpvt, double *tau, double *work, int *info);
 
int dgeqr2_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *info);
 
int dgeqrf_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *lwork, int *info);
 
int dgerfs_(char *trans, int *n, int *nrhs, 
	double *a, int *lda, double *af, int *ldaf, int *
	ipiv, double *b, int *ldb, double *x, int *ldx, 
	double *ferr, double *berr, double *work, int *iwork, 
	int *info);
 
int dgerq2_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *info);
 
int dgerqf_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *lwork, int *info);
 
int dgesc2_(int *n, double *a, int *lda, 
	double *rhs, int *ipiv, int *jpiv, double *scale);
 
int dgesdd_(char *jobz, int *m, int *n, double *
	a, int *lda, double *s, double *u, int *ldu, 
	double *vt, int *ldvt, double *work, int *lwork, 
	int *iwork, int *info);
 
int dgesv_(int *n, int *nrhs, double *a, int 
	*lda, int *ipiv, double *b, int *ldb, int *info);
 
int dgesvd_(char *jobu, char *jobvt, int *m, int *n, 
	double *a, int *lda, double *s, double *u, int *
	ldu, double *vt, int *ldvt, double *work, int *lwork, 
	int *info);
 
int dgesvx_(char *fact, char *trans, int *n, int *
	nrhs, double *a, int *lda, double *af, int *ldaf, 
	int *ipiv, char *equed, double *r__, double *c__, 
	double *b, int *ldb, double *x, int *ldx, double *
	rcond, double *ferr, double *berr, double *work, int *
	iwork, int *info);
 
int dgetc2_(int *n, double *a, int *lda, int 
	*ipiv, int *jpiv, int *info);
 
int dgetf2_(int *m, int *n, double *a, int *
	lda, int *ipiv, int *info);
 
int dgetrf_(int *m, int *n, double *a, int *
	lda, int *ipiv, int *info);
 
int dgetri_(int *n, double *a, int *lda, int 
	*ipiv, double *work, int *lwork, int *info);
 
int dgetrs_(char *trans, int *n, int *nrhs, 
	double *a, int *lda, int *ipiv, double *b, int *
	ldb, int *info);
 
int dggbak_(char *job, char *side, int *n, int *ilo, 
	int *ihi, double *lscale, double *rscale, int *m, 
	double *v, int *ldv, int *info);
 
int dggbal_(char *job, int *n, double *a, int *
	lda, double *b, int *ldb, int *ilo, int *ihi, 
	double *lscale, double *rscale, double *work, int *
	info);
 
int dgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	delctg, int *n, double *a, int *lda, double *b, 
	int *ldb, int *sdim, double *alphar, double *alphai, 
	double *beta, double *vsl, int *ldvsl, double *vsr, 
	int *ldvsr, double *work, int *lwork, int *bwork, 
	int *info);
 
int dggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	delctg, char *sense, int *n, double *a, int *lda, 
	double *b, int *ldb, int *sdim, double *alphar, 
	double *alphai, double *beta, double *vsl, int *ldvsl,
	 double *vsr, int *ldvsr, double *rconde, double *
	rcondv, double *work, int *lwork, int *iwork, int *
	liwork, int *bwork, int *info);
 
int dggev_(char *jobvl, char *jobvr, int *n, double *
	a, int *lda, double *b, int *ldb, double *alphar, 
	double *alphai, double *beta, double *vl, int *ldvl, 
	double *vr, int *ldvr, double *work, int *lwork, 
	int *info);
 
int dggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, int *n, double *a, int *lda, double *b, 
	int *ldb, double *alphar, double *alphai, double *
	beta, double *vl, int *ldvl, double *vr, int *ldvr, 
	int *ilo, int *ihi, double *lscale, double *rscale, 
	double *abnrm, double *bbnrm, double *rconde, double *
	rcondv, double *work, int *lwork, int *iwork, int *
	bwork, int *info);
 
int dggglm_(int *n, int *m, int *p, double *
	a, int *lda, double *b, int *ldb, double *d__, 
	double *x, double *y, double *work, int *lwork, 
	int *info);
 
int dgghrd_(char *compq, char *compz, int *n, int *
	ilo, int *ihi, double *a, int *lda, double *b, 
	int *ldb, double *q, int *ldq, double *z__, int *
	ldz, int *info);
 
int dgglse_(int *m, int *n, int *p, double *
	a, int *lda, double *b, int *ldb, double *c__, 
	double *d__, double *x, double *work, int *lwork, 
	int *info);
 
int dggqrf_(int *n, int *m, int *p, double *
	a, int *lda, double *taua, double *b, int *ldb, 
	double *taub, double *work, int *lwork, int *info);
 
int dggrqf_(int *m, int *p, int *n, double *
	a, int *lda, double *taua, double *b, int *ldb, 
	double *taub, double *work, int *lwork, int *info);
 
int dggsvd_(char *jobu, char *jobv, char *jobq, int *m, 
	int *n, int *p, int *k, int *l, double *a, 
	int *lda, double *b, int *ldb, double *alpha, 
	double *beta, double *u, int *ldu, double *v, int 
	*ldv, double *q, int *ldq, double *work, int *iwork, 
	int *info);
 
int dggsvp_(char *jobu, char *jobv, char *jobq, int *m, 
	int *p, int *n, double *a, int *lda, double *b, 
	int *ldb, double *tola, double *tolb, int *k, int 
	*l, double *u, int *ldu, double *v, int *ldv, 
	double *q, int *ldq, int *iwork, double *tau, 
	double *work, int *info);
 
int dgtcon_(char *norm, int *n, double *dl, 
	double *d__, double *du, double *du2, int *ipiv, 
	double *anorm, double *rcond, double *work, int *
	iwork, int *info);
 
int dgtrfs_(char *trans, int *n, int *nrhs, 
	double *dl, double *d__, double *du, double *dlf, 
	double *df, double *duf, double *du2, int *ipiv, 
	double *b, int *ldb, double *x, int *ldx, double *
	ferr, double *berr, double *work, int *iwork, int *
	info);
 
int dgtsv_(int *n, int *nrhs, double *dl, 
	double *d__, double *du, double *b, int *ldb, int 
	*info);
 
int dgtsvx_(char *fact, char *trans, int *n, int *
	nrhs, double *dl, double *d__, double *du, double *
	dlf, double *df, double *duf, double *du2, int *ipiv, 
	double *b, int *ldb, double *x, int *ldx, double *
	rcond, double *ferr, double *berr, double *work, int *
	iwork, int *info);
 
int dgttrf_(int *n, double *dl, double *d__, 
	double *du, double *du2, int *ipiv, int *info);
 
int dgttrs_(char *trans, int *n, int *nrhs, 
	double *dl, double *d__, double *du, double *du2, 
	int *ipiv, double *b, int *ldb, int *info);
 
int dgtts2_(int *itrans, int *n, int *nrhs, 
	double *dl, double *d__, double *du, double *du2, 
	int *ipiv, double *b, int *ldb);
 
int dhgeqz_(char *job, char *compq, char *compz, int *n, 
	int *ilo, int *ihi, double *a, int *lda, double *
	b, int *ldb, double *alphar, double *alphai, double *
	beta, double *q, int *ldq, double *z__, int *ldz, 
	double *work, int *lwork, int *info);
 
int dhsein_(char *side, char *eigsrc, char *initv, int *
	select, int *n, double *h__, int *ldh, double *wr, 
	double *wi, double *vl, int *ldvl, double *vr, 
	int *ldvr, int *mm, int *m, double *work, int *
	ifaill, int *ifailr, int *info);
 
int dhseqr_(char *job, char *compz, int *n, int *ilo,
	 int *ihi, double *h__, int *ldh, double *wr, 
	double *wi, double *z__, int *ldz, double *work, 
	int *lwork, int *info);
 
int dlabad_(double *small, double *large);
 
int dlabrd_(int *m, int *n, int *nb, double *
	a, int *lda, double *d__, double *e, double *tauq, 
	double *taup, double *x, int *ldx, double *y, int 
	*ldy);
 
int dlacon_(int *n, double *v, double *x, 
	int *isgn, double *est, int *kase);
 
int dlacpy_(char *uplo, int *m, int *n, double *
	a, int *lda, double *b, int *ldb);
 
int dladiv_(double *a, double *b, double *c__, 
	double *d__, double *p, double *q);
 
int dlae2_(double *a, double *b, double *c__, 
	double *rt1, double *rt2);
 
int dlaebz_(int *ijob, int *nitmax, int *n, 
	int *mmax, int *minp, int *nbmin, double *abstol, 
	double *reltol, double *pivmin, double *d__, double *
	e, double *e2, int *nval, double *ab, double *c__, 
	int *mout, int *nab, double *work, int *iwork, 
	int *info);
 
int dlaed0_(int *icompq, int *qsiz, int *n, 
	double *d__, double *e, double *q, int *ldq, 
	double *qstore, int *ldqs, double *work, int *iwork, 
	int *info);
 
int dlaed1_(int *n, double *d__, double *q, 
	int *ldq, int *indxq, double *rho, int *cutpnt, 
	double *work, int *iwork, int *info);
 
int dlaed2_(int *k, int *n, int *n1, double *
	d__, double *q, int *ldq, int *indxq, double *rho, 
	double *z__, double *dlamda, double *w, double *q2, 
	int *indx, int *indxc, int *indxp, int *coltyp, 
	int *info);
 
int dlaed3_(int *k, int *n, int *n1, double *
	d__, double *q, int *ldq, double *rho, double *dlamda,
	 double *q2, int *indx, int *ctot, double *w, 
	double *s, int *info);
 
int dlaed4_(int *n, int *i__, double *d__, 
	double *z__, double *delta, double *rho, double *dlam,
	 int *info);
 
int dlaed5_(int *i__, double *d__, double *z__, 
	double *delta, double *rho, double *dlam);
 
int dlaed6_(int *kniter, int *orgati, double *
	rho, double *d__, double *z__, double *finit, double *
	tau, int *info);
 
int dlaed7_(int *icompq, int *n, int *qsiz, 
	int *tlvls, int *curlvl, int *curpbm, double *d__, 
	double *q, int *ldq, int *indxq, double *rho, int 
	*cutpnt, double *qstore, int *qptr, int *prmptr, int *
	perm, int *givptr, int *givcol, double *givnum, 
	double *work, int *iwork, int *info);
 
int dlaed8_(int *icompq, int *k, int *n, int 
	*qsiz, double *d__, double *q, int *ldq, int *indxq, 
	double *rho, int *cutpnt, double *z__, double *dlamda,
	 double *q2, int *ldq2, double *w, int *perm, int 
	*givptr, int *givcol, double *givnum, int *indxp, int 
	*indx, int *info);
 
int dlaed9_(int *k, int *kstart, int *kstop, 
	int *n, double *d__, double *q, int *ldq, double *
	rho, double *dlamda, double *w, double *s, int *lds, 
	int *info);
 
int dlaeda_(int *n, int *tlvls, int *curlvl, 
	int *curpbm, int *prmptr, int *perm, int *givptr, 
	int *givcol, double *givnum, double *q, int *qptr, 
	double *z__, double *ztemp, int *info);
 
int dlaein_(int *rightv, int *noinit, int *n, 
	double *h__, int *ldh, double *wr, double *wi, 
	double *vr, double *vi, double *b, int *ldb, 
	double *work, double *eps3, double *smlnum, double *
	bignum, int *info);
 
int dlaev2_(double *a, double *b, double *c__, 
	double *rt1, double *rt2, double *cs1, double *sn1);
 
int dlaexc_(int *wantq, int *n, double *t, 
	int *ldt, double *q, int *ldq, int *j1, int *n1, 
	int *n2, double *work, int *info);
 
int dlag2_(double *a, int *lda, double *b, 
	int *ldb, double *safmin, double *scale1, double *
	scale2, double *wr1, double *wr2, double *wi);
 
int dlags2_(int *upper, double *a1, double *a2, 
	double *a3, double *b1, double *b2, double *b3, 
	double *csu, double *snu, double *csv, double *snv, 
	double *csq, double *snq);
 
int dlagtf_(int *n, double *a, double *lambda, 
	double *b, double *c__, double *tol, double *d__, 
	int *in, int *info);
 
int dlagtm_(char *trans, int *n, int *nrhs, 
	double *alpha, double *dl, double *d__, double *du, 
	double *x, int *ldx, double *beta, double *b, int 
	*ldb);
 
int dlagts_(int *job, int *n, double *a, 
	double *b, double *c__, double *d__, int *in, 
	double *y, double *tol, int *info);
 
int dlagv2_(double *a, int *lda, double *b, 
	int *ldb, double *alphar, double *alphai, double *
	beta, double *csl, double *snl, double *csr, double *
	snr);
 
int dlahqr_(int *wantt, int *wantz, int *n, 
	int *ilo, int *ihi, double *h__, int *ldh, double 
	*wr, double *wi, int *iloz, int *ihiz, double *z__, 
	int *ldz, int *info);
 
int dlahrd_(int *n, int *k, int *nb, double *
	a, int *lda, double *tau, double *t, int *ldt, 
	double *y, int *ldy);
 
int dlaic1_(int *job, int *j, double *x, 
	double *sest, double *w, double *gamma, double *
	sestpr, double *s, double *c__);
 
int dlaln2_(int *ltrans, int *na, int *nw, 
	double *smin, double *ca, double *a, int *lda, 
	double *d1, double *d2, double *b, int *ldb, 
	double *wr, double *wi, double *x, int *ldx, 
	double *scale, double *xnorm, int *info);
 
int dlals0_(int *icompq, int *nl, int *nr, 
	int *sqre, int *nrhs, double *b, int *ldb, double 
	*bx, int *ldbx, int *perm, int *givptr, int *givcol, 
	int *ldgcol, double *givnum, int *ldgnum, double *
	poles, double *difl, double *difr, double *z__, int *
	k, double *c__, double *s, double *work, int *info);
 
int dlalsa_(int *icompq, int *smlsiz, int *n, 
	int *nrhs, double *b, int *ldb, double *bx, int *
	ldbx, double *u, int *ldu, double *vt, int *k, 
	double *difl, double *difr, double *z__, double *
	poles, int *givptr, int *givcol, int *ldgcol, int *
	perm, double *givnum, double *c__, double *s, double *
	work, int *iwork, int *info);
 
int dlalsd_(char *uplo, int *smlsiz, int *n, int 
	*nrhs, double *d__, double *e, double *b, int *ldb, 
	double *rcond, int *rank, double *work, int *iwork, 
	int *info);
 
int dlamc1_(int *beta, int *t, int *rnd, int 
	*ieee1);
 
int dlamc2_(int *beta, int *t, int *rnd, 
	double *eps, int *emin, double *rmin, int *emax, 
	double *rmax);
 
int dlamc4_(int *emin, double *start, int *base);
 
int dlamc5_(int *beta, int *p, int *emin, 
	int *ieee, int *emax, double *rmax);
 
int dlamrg_(int *n1, int *n2, double *a, int 
	*dtrd1, int *dtrd2, int *index);
 
int dlanv2_(double *a, double *b, double *c__, 
	double *d__, double *rt1r, double *rt1i, double *rt2r,
	 double *rt2i, double *cs, double *sn);
 
int dlapll_(int *n, double *x, int *incx, 
	double *y, int *incy, double *ssmin);
 
int dlapmt_(int *forwrd, int *m, int *n, 
	double *x, int *ldx, int *k);
 
int dlaqgb_(int *m, int *n, int *kl, int *ku,
	 double *ab, int *ldab, double *r__, double *c__, 
	double *rowcnd, double *colcnd, double *amax, char *equed);
 
int dlaqge_(int *m, int *n, double *a, int *
	lda, double *r__, double *c__, double *rowcnd, double 
	*colcnd, double *amax, char *equed);
 
int dlaqp2_(int *m, int *n, int *offset, 
	double *a, int *lda, int *jpvt, double *tau, 
	double *vn1, double *vn2, double *work);
 
int dlaqps_(int *m, int *n, int *offset, int 
	*nb, int *kb, double *a, int *lda, int *jpvt, 
	double *tau, double *vn1, double *vn2, double *auxv, 
	double *f, int *ldf);
 
int dlaqsb_(char *uplo, int *n, int *kd, double *
	ab, int *ldab, double *s, double *scond, double *amax,
	 char *equed);
 
int dlaqsp_(char *uplo, int *n, double *ap, 
	double *s, double *scond, double *amax, char *equed);
 
int dlaqsy_(char *uplo, int *n, double *a, int *
	lda, double *s, double *scond, double *amax, char *equed);
 
int dlaqtr_(int *ltran, int *lfloat, int *n, 
	double *t, int *ldt, double *b, double *w, double 
	*scale, double *x, double *work, int *info);
 
int dlar1v_(int *n, int *b1, int *bn, double 
	*sigma, double *d__, double *l, double *ld, double *
	lld, double *gersch, double *z__, double *ztz, double 
	*mingma, int *r__, int *isuppz, double *work);
 
int dlar2v_(int *n, double *x, double *y, 
	double *z__, int *incx, double *c__, double *s, 
	int *incc);
 
int dlarf_(char *side, int *m, int *n, double *v,
	 int *incv, double *tau, double *c__, int *ldc, 
	double *work);
 
int dlarfb_(char *side, char *trans, char *direct, char *
	storev, int *m, int *n, int *k, double *v, int *
	ldv, double *t, int *ldt, double *c__, int *ldc, 
	double *work, int *ldwork);
 
int dlarfg_(int *n, double *alpha, double *x, 
	int *incx, double *tau);
 
int dlarft_(char *direct, char *storev, int *n, int *
	k, double *v, int *ldv, double *tau, double *t, 
	int *ldt);
 
int dlarfx_(char *side, int *m, int *n, double *
	v, double *tau, double *c__, int *ldc, double *work);
 
int dlargv_(int *n, double *x, int *incx, 
	double *y, int *incy, double *c__, int *incc);
 
int dlarnv_(int *idist, int *iseed, int *n, 
	double *x);
 
int dlarrb_(int *n, double *d__, double *l, 
	double *ld, double *lld, int *ifirst, int *ilast, 
	double *sigma, double *reltol, double *w, double *
	wgap, double *werr, double *work, int *iwork, int *
	info);
 
int dlarre_(int *n, double *d__, double *e, 
	double *tol, int *nsplit, int *isplit, int *m, 
	double *w, double *woff, double *gersch, double *work,
	 int *info);
 
int dlarrf_(int *n, double *d__, double *l, 
	double *ld, double *lld, int *ifirst, int *ilast, 
	double *w, double *dplus, double *lplus, double *work,
	 int *iwork, int *info);
 
int dlarrv_(int *n, double *d__, double *l, 
	int *isplit, int *m, double *w, int *iblock, 
	double *gersch, double *tol, double *z__, int *ldz, 
	int *isuppz, double *work, int *iwork, int *info);
 
int dlartg_(double *f, double *g, double *cs, 
	double *sn, double *r__);
 
int dlartv_(int *n, double *x, int *incx, 
	double *y, int *incy, double *c__, double *s, int 
	*incc);
 
int dlaruv_(int *iseed, int *n, double *x);
 
int dlarz_(char *side, int *m, int *n, int *l, 
	double *v, int *incv, double *tau, double *c__, 
	int *ldc, double *work);
 
int dlarzb_(char *side, char *trans, char *direct, char *
	storev, int *m, int *n, int *k, int *l, double *v,
	 int *ldv, double *t, int *ldt, double *c__, int *
	ldc, double *work, int *ldwork);
 
int dlarzt_(char *direct, char *storev, int *n, int *
	k, double *v, int *ldv, double *tau, double *t, 
	int *ldt);
 
int dlas2_(double *f, double *g, double *h__, 
	double *ssmin, double *ssmax);
 
int dlascl_(char *type__, int *kl, int *ku, 
	double *cfrom, double *cto, int *m, int *n, 
	double *a, int *lda, int *info);
 
int dlasd0_(int *n, int *sqre, double *d__, 
	double *e, double *u, int *ldu, double *vt, int *
	ldvt, int *smlsiz, int *iwork, double *work, int *
	info);
 
int dlasd1_(int *nl, int *nr, int *sqre, 
	double *d__, double *alpha, double *beta, double *u, 
	int *ldu, double *vt, int *ldvt, int *idxq, int *
	iwork, double *work, int *info);
 
int dlasd2_(int *nl, int *nr, int *sqre, int 
	*k, double *d__, double *z__, double *alpha, double *
	beta, double *u, int *ldu, double *vt, int *ldvt, 
	double *dsigma, double *u2, int *ldu2, double *vt2, 
	int *ldvt2, int *idxp, int *idx, int *idxc, int *
	idxq, int *coltyp, int *info);
 
int dlasd3_(int *nl, int *nr, int *sqre, int 
	*k, double *d__, double *q, int *ldq, double *dsigma, 
	double *u, int *ldu, double *u2, int *ldu2, 
	double *vt, int *ldvt, double *vt2, int *ldvt2, 
	int *idxc, int *ctot, double *z__, int *info);
 
int dlasd4_(int *n, int *i__, double *d__, 
	double *z__, double *delta, double *rho, double *
	sigma, double *work, int *info);
 
int dlasd5_(int *i__, double *d__, double *z__, 
	double *delta, double *rho, double *dsigma, double *
	work);
 
int dlasd6_(int *icompq, int *nl, int *nr, 
	int *sqre, double *d__, double *vf, double *vl, 
	double *alpha, double *beta, int *idxq, int *perm, 
	int *givptr, int *givcol, int *ldgcol, double *givnum,
	 int *ldgnum, double *poles, double *difl, double *
	difr, double *z__, int *k, double *c__, double *s, 
	double *work, int *iwork, int *info);
 
int dlasd7_(int *icompq, int *nl, int *nr, 
	int *sqre, int *k, double *d__, double *z__, 
	double *zw, double *vf, double *vfw, double *vl, 
	double *vlw, double *alpha, double *beta, double *
	dsigma, int *idx, int *idxp, int *idxq, int *perm, 
	int *givptr, int *givcol, int *ldgcol, double *givnum,
	 int *ldgnum, double *c__, double *s, int *info);
 
int dlasd8_(int *icompq, int *k, double *d__, 
	double *z__, double *vf, double *vl, double *difl, 
	double *difr, int *lddifr, double *dsigma, double *
	work, int *info);
 
int dlasd9_(int *icompq, int *ldu, int *k, 
	double *d__, double *z__, double *vf, double *vl, 
	double *difl, double *difr, double *dsigma, double *
	work, int *info);
 
int dlasda_(int *icompq, int *smlsiz, int *n, 
	int *sqre, double *d__, double *e, double *u, int 
	*ldu, double *vt, int *k, double *difl, double *difr, 
	double *z__, double *poles, int *givptr, int *givcol, 
	int *ldgcol, int *perm, double *givnum, double *c__, 
	double *s, double *work, int *iwork, int *info);
 
int dlasdq_(char *uplo, int *sqre, int *n, int *
	ncvt, int *nru, int *ncc, double *d__, double *e, 
	double *vt, int *ldvt, double *u, int *ldu, 
	double *c__, int *ldc, double *work, int *info);
 
int dlasdt_(int *n, int *lvl, int *nd, int *
	inode, int *ndiml, int *ndimr, int *msub);
 
int dlaset_(char *uplo, int *m, int *n, double *
	alpha, double *beta, double *a, int *lda);
 
int dlasq1_(int *n, double *d__, double *e, 
	double *work, int *info);
 
int dlasq2_(int *n, double *z__, int *info);
 
int dlasq3_(int *i0, int *n0, double *z__, 
	int *pp, double *dmin__, double *sigma, double *desig,
	 double *qmax, int *nfail, int *iter, int *ndiv, 
	int *ieee);
 
int dlasq4_(int *i0, int *n0, double *z__, 
	int *pp, int *n0in, double *dmin__, double *dmin1, 
	double *dmin2, double *dn, double *dn1, double *dn2, 
	double *tau, int *ttype);
 
int dlasq5_(int *i0, int *n0, double *z__, 
	int *pp, double *tau, double *dmin__, double *dmin1, 
	double *dmin2, double *dn, double *dnm1, double *dnm2,
	 int *ieee);
 
int dlasq6_(int *i0, int *n0, double *z__, 
	int *pp, double *dmin__, double *dmin1, double *dmin2,
	 double *dn, double *dnm1, double *dnm2);
 
int dlasr_(char *side, char *pivot, char *direct, int *m,
	 int *n, double *c__, double *s, double *a, int *
	lda);
 
int dlasrt_(char *id, int *n, double *d__, int *
	info);
 
int dlassq_(int *n, double *x, int *incx, 
	double *scale, double *sumsq);
 
int dlasv2_(double *f, double *g, double *h__, 
	double *ssmin, double *ssmax, double *snr, double *
	csr, double *snl, double *csl);
 
int dlaswp_(int *n, double *a, int *lda, int 
	*k1, int *k2, int *ipiv, int *incx);
 
int dlasy2_(int *ltranl, int *ltranr, int *isgn, 
	int *n1, int *n2, double *tl, int *ldtl, double *
	tr, int *ldtr, double *b, int *ldb, double *scale, 
	double *x, int *ldx, double *xnorm, int *info);
 
int dlasyf_(char *uplo, int *n, int *nb, int *kb,
	 double *a, int *lda, int *ipiv, double *w, int *
	ldw, int *info);
 
int dlatbs_(char *uplo, char *trans, char *diag, char *
	normin, int *n, int *kd, double *ab, int *ldab, 
	double *x, double *scale, double *cnorm, int *info);
 
int dlatdf_(int *ijob, int *n, double *z__, 
	int *ldz, double *rhs, double *rdsum, double *rdscal, 
	int *ipiv, int *jpiv);
 
int dlatps_(char *uplo, char *trans, char *diag, char *
	normin, int *n, double *ap, double *x, double *scale, 
	double *cnorm, int *info);
 
int dlatrd_(char *uplo, int *n, int *nb, double *
	a, int *lda, double *e, double *tau, double *w, 
	int *ldw);
 
int dlatrs_(char *uplo, char *trans, char *diag, char *
	normin, int *n, double *a, int *lda, double *x, 
	double *scale, double *cnorm, int *info);
 
int dlatrz_(int *m, int *n, int *l, double *
	a, int *lda, double *tau, double *work);
 
int dlatzm_(char *side, int *m, int *n, double *
	v, int *incv, double *tau, double *c1, double *c2, 
	int *ldc, double *work);
 
int dlauu2_(char *uplo, int *n, double *a, int *
	lda, int *info);
 
int dlauum_(char *uplo, int *n, double *a, int *
	lda, int *info);
 
int dopgtr_(char *uplo, int *n, double *ap, 
	double *tau, double *q, int *ldq, double *work, 
	int *info);
 
int dopmtr_(char *side, char *uplo, char *trans, int *m, 
	int *n, double *ap, double *tau, double *c__, int 
	*ldc, double *work, int *info);
 
int dorg2l_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *info);
 
int dorg2r_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *info);
 
int dorgbr_(char *vect, int *m, int *n, int *k, 
	double *a, int *lda, double *tau, double *work, 
	int *lwork, int *info);
 
int dorghr_(int *n, int *ilo, int *ihi, 
	double *a, int *lda, double *tau, double *work, 
	int *lwork, int *info);
 
int dorgl2_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *info);
 
int dorglq_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *lwork, 
	int *info);
 
int dorgql_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *lwork, 
	int *info);
 
int dorgqr_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *lwork, 
	int *info);
 
int dorgr2_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *info);
 
int dorgrq_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *lwork, 
	int *info);
 
int dorgtr_(char *uplo, int *n, double *a, int *
	lda, double *tau, double *work, int *lwork, int *info);
 
int dorm2l_(char *side, char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *info);
 
int dorm2r_(char *side, char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *info);
 
int dormbr_(char *vect, char *side, char *trans, int *m, 
	int *n, int *k, double *a, int *lda, double *tau, 
	double *c__, int *ldc, double *work, int *lwork, 
	int *info);
 
int dormhr_(char *side, char *trans, int *m, int *n, 
	int *ilo, int *ihi, double *a, int *lda, double *
	tau, double *c__, int *ldc, double *work, int *lwork, 
	int *info);
 
int dorml2_(char *side, char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *info);
 
int dormlq_(char *side, char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *lwork, int *info);
 
int dormql_(char *side, char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *lwork, int *info);
 
int dormqr_(char *side, char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *lwork, int *info);
 
int dormr2_(char *side, char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *info);
 
int dormr3_(char *side, char *trans, int *m, int *n, 
	int *k, int *l, double *a, int *lda, double *tau, 
	double *c__, int *ldc, double *work, int *info);
 
int dormrq_(char *side, char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *lwork, int *info);
 
int dormrz_(char *side, char *trans, int *m, int *n, 
	int *k, int *l, double *a, int *lda, double *tau, 
	double *c__, int *ldc, double *work, int *lwork, 
	int *info);
 
int dormtr_(char *side, char *uplo, char *trans, int *m, 
	int *n, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *lwork, int *info);
 
int dpbcon_(char *uplo, int *n, int *kd, double *
	ab, int *ldab, double *anorm, double *rcond, double *
	work, int *iwork, int *info);
 
int dpbequ_(char *uplo, int *n, int *kd, double *
	ab, int *ldab, double *s, double *scond, double *amax,
	 int *info);
 
int dpbrfs_(char *uplo, int *n, int *kd, int *
	nrhs, double *ab, int *ldab, double *afb, int *ldafb, 
	double *b, int *ldb, double *x, int *ldx, double *
	ferr, double *berr, double *work, int *iwork, int *
	info);
 
int dpbstf_(char *uplo, int *n, int *kd, double *
	ab, int *ldab, int *info);
 
int dpbsv_(char *uplo, int *n, int *kd, int *
	nrhs, double *ab, int *ldab, double *b, int *ldb, 
	int *info);
 
int dpbsvx_(char *fact, char *uplo, int *n, int *kd, 
	int *nrhs, double *ab, int *ldab, double *afb, 
	int *ldafb, char *equed, double *s, double *b, int *
	ldb, double *x, int *ldx, double *rcond, double *ferr,
	 double *berr, double *work, int *iwork, int *info);
 
int dpbtf2_(char *uplo, int *n, int *kd, double *
	ab, int *ldab, int *info);
 
int dpbtrf_(char *uplo, int *n, int *kd, double *
	ab, int *ldab, int *info);
 
int dpbtrs_(char *uplo, int *n, int *kd, int *
	nrhs, double *ab, int *ldab, double *b, int *ldb, 
	int *info);
 
int dpocon_(char *uplo, int *n, double *a, int *
	lda, double *anorm, double *rcond, double *work, int *
	iwork, int *info);
 
int dpoequ_(int *n, double *a, int *lda, 
	double *s, double *scond, double *amax, int *info);
 
int dporfs_(char *uplo, int *n, int *nrhs, 
	double *a, int *lda, double *af, int *ldaf, 
	double *b, int *ldb, double *x, int *ldx, double *
	ferr, double *berr, double *work, int *iwork, int *
	info);
 
int dposv_(char *uplo, int *n, int *nrhs, double 
	*a, int *lda, double *b, int *ldb, int *info);
 
int dposvx_(char *fact, char *uplo, int *n, int *
	nrhs, double *a, int *lda, double *af, int *ldaf, 
	char *equed, double *s, double *b, int *ldb, double *
	x, int *ldx, double *rcond, double *ferr, double *
	berr, double *work, int *iwork, int *info);
 
int dpotf2_(char *uplo, int *n, double *a, int *
	lda, int *info);
 
int dpotrf_(char *uplo, int *n, double *a, int *
	lda, int *info);
 
int dpotri_(char *uplo, int *n, double *a, int *
	lda, int *info);
 
int dpotrs_(char *uplo, int *n, int *nrhs, 
	double *a, int *lda, double *b, int *ldb, int *
	info);
 
int dppcon_(char *uplo, int *n, double *ap, 
	double *anorm, double *rcond, double *work, int *
	iwork, int *info);
 
int dppequ_(char *uplo, int *n, double *ap, 
	double *s, double *scond, double *amax, int *info);
 
int dpprfs_(char *uplo, int *n, int *nrhs, 
	double *ap, double *afp, double *b, int *ldb, 
	double *x, int *ldx, double *ferr, double *berr, 
	double *work, int *iwork, int *info);
 
int dppsv_(char *uplo, int *n, int *nrhs, double 
	*ap, double *b, int *ldb, int *info);
 
int dppsvx_(char *fact, char *uplo, int *n, int *
	nrhs, double *ap, double *afp, char *equed, double *s, 
	double *b, int *ldb, double *x, int *ldx, double *
	rcond, double *ferr, double *berr, double *work, int *
	iwork, int *info);
 
int dpptrf_(char *uplo, int *n, double *ap, int *
	info);
 
int dpptri_(char *uplo, int *n, double *ap, int *
	info);
 
int dpptrs_(char *uplo, int *n, int *nrhs, 
	double *ap, double *b, int *ldb, int *info);
 
int dptcon_(int *n, double *d__, double *e, 
	double *anorm, double *rcond, double *work, int *info);
 
int dpteqr_(char *compz, int *n, double *d__, 
	double *e, double *z__, int *ldz, double *work, 
	int *info);
 
int dptrfs_(int *n, int *nrhs, double *d__, 
	double *e, double *df, double *ef, double *b, int 
	*ldb, double *x, int *ldx, double *ferr, double *berr,
	 double *work, int *info);
 
int dptsv_(int *n, int *nrhs, double *d__, 
	double *e, double *b, int *ldb, int *info);
 
int dptsvx_(char *fact, int *n, int *nrhs, 
	double *d__, double *e, double *df, double *ef, 
	double *b, int *ldb, double *x, int *ldx, double *
	rcond, double *ferr, double *berr, double *work, int *
	info);
 
int dpttrf_(int *n, double *d__, double *e, 
	int *info);
 
int dpttrs_(int *n, int *nrhs, double *d__, 
	double *e, double *b, int *ldb, int *info);
 
int dptts2_(int *n, int *nrhs, double *d__, 
	double *e, double *b, int *ldb);
 
int drscl_(int *n, double *sa, double *sx, 
	int *incx);
 
int dsbev_(char *jobz, char *uplo, int *n, int *kd, 
	double *ab, int *ldab, double *w, double *z__, 
	int *ldz, double *work, int *info);
 
int dsbevd_(char *jobz, char *uplo, int *n, int *kd, 
	double *ab, int *ldab, double *w, double *z__, 
	int *ldz, double *work, int *lwork, int *iwork, 
	int *liwork, int *info);
 
int dsbevx_(char *jobz, char *range, char *uplo, int *n, 
	int *kd, double *ab, int *ldab, double *q, int *
	ldq, double *vl, double *vu, int *il, int *iu, 
	double *abstol, int *m, double *w, double *z__, 
	int *ldz, double *work, int *iwork, int *ifail, 
	int *info);
 
int dsbgst_(char *vect, char *uplo, int *n, int *ka, 
	int *kb, double *ab, int *ldab, double *bb, int *
	ldbb, double *x, int *ldx, double *work, int *info);
 
int dsbgv_(char *jobz, char *uplo, int *n, int *ka, 
	int *kb, double *ab, int *ldab, double *bb, int *
	ldbb, double *w, double *z__, int *ldz, double *work, 
	int *info);
 
int dsbgvd_(char *jobz, char *uplo, int *n, int *ka, 
	int *kb, double *ab, int *ldab, double *bb, int *
	ldbb, double *w, double *z__, int *ldz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
int dsbgvx_(char *jobz, char *range, char *uplo, int *n, 
	int *ka, int *kb, double *ab, int *ldab, double *
	bb, int *ldbb, double *q, int *ldq, double *vl, 
	double *vu, int *il, int *iu, double *abstol, int 
	*m, double *w, double *z__, int *ldz, double *work, 
	int *iwork, int *ifail, int *info);
 
int dsbtrd_(char *vect, char *uplo, int *n, int *kd, 
	double *ab, int *ldab, double *d__, double *e, 
	double *q, int *ldq, double *work, int *info);
 
int dspcon_(char *uplo, int *n, double *ap, int *
	ipiv, double *anorm, double *rcond, double *work, int 
	*iwork, int *info);
 
int dspev_(char *jobz, char *uplo, int *n, double *
	ap, double *w, double *z__, int *ldz, double *work, 
	int *info);
 
int dspevd_(char *jobz, char *uplo, int *n, double *
	ap, double *w, double *z__, int *ldz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
int dspevx_(char *jobz, char *range, char *uplo, int *n, 
	double *ap, double *vl, double *vu, int *il, int *
	iu, double *abstol, int *m, double *w, double *z__, 
	int *ldz, double *work, int *iwork, int *ifail, 
	int *info);
 
int dspgst_(int *itype, char *uplo, int *n, 
	double *ap, double *bp, int *info);
 
int dspgv_(int *itype, char *jobz, char *uplo, int *
	n, double *ap, double *bp, double *w, double *z__, 
	int *ldz, double *work, int *info);
 
int dspgvd_(int *itype, char *jobz, char *uplo, int *
	n, double *ap, double *bp, double *w, double *z__, 
	int *ldz, double *work, int *lwork, int *iwork, 
	int *liwork, int *info);
 
int dspgvx_(int *itype, char *jobz, char *range, char *
	uplo, int *n, double *ap, double *bp, double *vl, 
	double *vu, int *il, int *iu, double *abstol, int 
	*m, double *w, double *z__, int *ldz, double *work, 
	int *iwork, int *ifail, int *info);
 
int dsprfs_(char *uplo, int *n, int *nrhs, 
	double *ap, double *afp, int *ipiv, double *b, 
	int *ldb, double *x, int *ldx, double *ferr, 
	double *berr, double *work, int *iwork, int *info);
 
int dspsv_(char *uplo, int *n, int *nrhs, double 
	*ap, int *ipiv, double *b, int *ldb, int *info);
 
int dspsvx_(char *fact, char *uplo, int *n, int *
	nrhs, double *ap, double *afp, int *ipiv, double *b, 
	int *ldb, double *x, int *ldx, double *rcond, 
	double *ferr, double *berr, double *work, int *iwork, 
	int *info);
 
int dsptrd_(char *uplo, int *n, double *ap, 
	double *d__, double *e, double *tau, int *info);
 
int dsptrf_(char *uplo, int *n, double *ap, int *
	ipiv, int *info);
 
int dsptri_(char *uplo, int *n, double *ap, int *
	ipiv, double *work, int *info);
 
int dsptrs_(char *uplo, int *n, int *nrhs, 
	double *ap, int *ipiv, double *b, int *ldb, int *
	info);
 
int dstebz_(char *range, char *order, int *n, double 
	*vl, double *vu, int *il, int *iu, double *abstol, 
	double *d__, double *e, int *m, int *nsplit, 
	double *w, int *iblock, int *isplit, double *work, 
	int *iwork, int *info);
 
int dstedc_(char *compz, int *n, double *d__, 
	double *e, double *z__, int *ldz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
int dstegr_(char *jobz, char *range, int *n, double *
	d__, double *e, double *vl, double *vu, int *il, 
	int *iu, double *abstol, int *m, double *w, 
	double *z__, int *ldz, int *isuppz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
int dstein_(int *n, double *d__, double *e, 
	int *m, double *w, int *iblock, int *isplit, 
	double *z__, int *ldz, double *work, int *iwork, 
	int *ifail, int *info);
 
int dsteqr_(char *compz, int *n, double *d__, 
	double *e, double *z__, int *ldz, double *work, 
	int *info);
 
int dsterf_(int *n, double *d__, double *e, 
	int *info);
 
int dstev_(char *jobz, int *n, double *d__, 
	double *e, double *z__, int *ldz, double *work, 
	int *info);
 
int dstevd_(char *jobz, int *n, double *d__, 
	double *e, double *z__, int *ldz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
int dstevr_(char *jobz, char *range, int *n, double *
	d__, double *e, double *vl, double *vu, int *il, 
	int *iu, double *abstol, int *m, double *w, 
	double *z__, int *ldz, int *isuppz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
int dstevx_(char *jobz, char *range, int *n, double *
	d__, double *e, double *vl, double *vu, int *il, 
	int *iu, double *abstol, int *m, double *w, 
	double *z__, int *ldz, double *work, int *iwork, 
	int *ifail, int *info);
 
int dsycon_(char *uplo, int *n, double *a, int *
	lda, int *ipiv, double *anorm, double *rcond, double *
	work, int *iwork, int *info);
 
int dsyev_(char *jobz, char *uplo, int *n, double *a,
	 int *lda, double *w, double *work, int *lwork, 
	int *info);
 
int dsyevd_(char *jobz, char *uplo, int *n, double *
	a, int *lda, double *w, double *work, int *lwork, 
	int *iwork, int *liwork, int *info);
 
int dsyevr_(char *jobz, char *range, char *uplo, int *n, 
	double *a, int *lda, double *vl, double *vu, int *
	il, int *iu, double *abstol, int *m, double *w, 
	double *z__, int *ldz, int *isuppz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
int dsyevx_(char *jobz, char *range, char *uplo, int *n, 
	double *a, int *lda, double *vl, double *vu, int *
	il, int *iu, double *abstol, int *m, double *w, 
	double *z__, int *ldz, double *work, int *lwork, 
	int *iwork, int *ifail, int *info);
 
int dsygs2_(int *itype, char *uplo, int *n, 
	double *a, int *lda, double *b, int *ldb, int *
	info);
 
int dsygst_(int *itype, char *uplo, int *n, 
	double *a, int *lda, double *b, int *ldb, int *
	info);
 
int dsygv_(int *itype, char *jobz, char *uplo, int *
	n, double *a, int *lda, double *b, int *ldb, 
	double *w, double *work, int *lwork, int *info);
 
int dsygvd_(int *itype, char *jobz, char *uplo, int *
	n, double *a, int *lda, double *b, int *ldb, 
	double *w, double *work, int *lwork, int *iwork, 
	int *liwork, int *info);
 
int dsygvx_(int *itype, char *jobz, char *range, char *
	uplo, int *n, double *a, int *lda, double *b, int 
	*ldb, double *vl, double *vu, int *il, int *iu, 
	double *abstol, int *m, double *w, double *z__, 
	int *ldz, double *work, int *lwork, int *iwork, 
	int *ifail, int *info);
 
int dsyrfs_(char *uplo, int *n, int *nrhs, 
	double *a, int *lda, double *af, int *ldaf, int *
	ipiv, double *b, int *ldb, double *x, int *ldx, 
	double *ferr, double *berr, double *work, int *iwork, 
	int *info);
 
int dsysv_(char *uplo, int *n, int *nrhs, double 
	*a, int *lda, int *ipiv, double *b, int *ldb, 
	double *work, int *lwork, int *info);
 
int dsysvx_(char *fact, char *uplo, int *n, int *
	nrhs, double *a, int *lda, double *af, int *ldaf, 
	int *ipiv, double *b, int *ldb, double *x, int *
	ldx, double *rcond, double *ferr, double *berr, 
	double *work, int *lwork, int *iwork, int *info);
 
int dsytd2_(char *uplo, int *n, double *a, int *
	lda, double *d__, double *e, double *tau, int *info);
 
int dsytf2_(char *uplo, int *n, double *a, int *
	lda, int *ipiv, int *info);
 
int dsytrd_(char *uplo, int *n, double *a, int *
	lda, double *d__, double *e, double *tau, double *
	work, int *lwork, int *info);
 
int dsytrf_(char *uplo, int *n, double *a, int *
	lda, int *ipiv, double *work, int *lwork, int *info);
 
int dsytri_(char *uplo, int *n, double *a, int *
	lda, int *ipiv, double *work, int *info);
 
int dsytrs_(char *uplo, int *n, int *nrhs, 
	double *a, int *lda, int *ipiv, double *b, int *
	ldb, int *info);
 
int dtbcon_(char *norm, char *uplo, char *diag, int *n, 
	int *kd, double *ab, int *ldab, double *rcond, 
	double *work, int *iwork, int *info);
 
int dtbrfs_(char *uplo, char *trans, char *diag, int *n, 
	int *kd, int *nrhs, double *ab, int *ldab, double 
	*b, int *ldb, double *x, int *ldx, double *ferr, 
	double *berr, double *work, int *iwork, int *info);
 
int dtbtrs_(char *uplo, char *trans, char *diag, int *n, 
	int *kd, int *nrhs, double *ab, int *ldab, double 
	*b, int *ldb, int *info);
 
int dtgevc_(char *side, char *howmny, int *select, 
	int *n, double *a, int *lda, double *b, int *ldb, 
	double *vl, int *ldvl, double *vr, int *ldvr, int 
	*mm, int *m, double *work, int *info);
 
int dtgex2_(int *wantq, int *wantz, int *n, 
	double *a, int *lda, double *b, int *ldb, double *
	q, int *ldq, double *z__, int *ldz, int *j1, int *
	n1, int *n2, double *work, int *lwork, int *info);
 
int dtgexc_(int *wantq, int *wantz, int *n, 
	double *a, int *lda, double *b, int *ldb, double *
	q, int *ldq, double *z__, int *ldz, int *ifst, 
	int *ilst, double *work, int *lwork, int *info);
 
int dtgsen_(int *ijob, int *wantq, int *wantz, 
	int *select, int *n, double *a, int *lda, double *
	b, int *ldb, double *alphar, double *alphai, double *
	beta, double *q, int *ldq, double *z__, int *ldz, 
	int *m, double *pl, double *pr, double *dif, 
	double *work, int *lwork, int *iwork, int *liwork, 
	int *info);
 
int dtgsja_(char *jobu, char *jobv, char *jobq, int *m, 
	int *p, int *n, int *k, int *l, double *a, 
	int *lda, double *b, int *ldb, double *tola, 
	double *tolb, double *alpha, double *beta, double *u, 
	int *ldu, double *v, int *ldv, double *q, int *
	ldq, double *work, int *ncycle, int *info);
 
int dtgsna_(char *job, char *howmny, int *select, 
	int *n, double *a, int *lda, double *b, int *ldb, 
	double *vl, int *ldvl, double *vr, int *ldvr, 
	double *s, double *dif, int *mm, int *m, double *
	work, int *lwork, int *iwork, int *info);
 
int dtgsy2_(char *trans, int *ijob, int *m, int *
	n, double *a, int *lda, double *b, int *ldb, 
	double *c__, int *ldc, double *d__, int *ldd, 
	double *e, int *lde, double *f, int *ldf, double *
	scale, double *rdsum, double *rdscal, int *iwork, int 
	*pq, int *info);
 
int dtgsyl_(char *trans, int *ijob, int *m, int *
	n, double *a, int *lda, double *b, int *ldb, 
	double *c__, int *ldc, double *d__, int *ldd, 
	double *e, int *lde, double *f, int *ldf, double *
	scale, double *dif, double *work, int *lwork, int *
	iwork, int *info);
 
int dtpcon_(char *norm, char *uplo, char *diag, int *n, 
	double *ap, double *rcond, double *work, int *iwork, 
	int *info);
 
int dtprfs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, double *ap, double *b, int *ldb, 
	double *x, int *ldx, double *ferr, double *berr, 
	double *work, int *iwork, int *info);
 
int dtptri_(char *uplo, char *diag, int *n, double *
	ap, int *info);
 
int dtptrs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, double *ap, double *b, int *ldb, int *
	info);
 
int dtrcon_(char *norm, char *uplo, char *diag, int *n, 
	double *a, int *lda, double *rcond, double *work, 
	int *iwork, int *info);
 
int dtrevc_(char *side, char *howmny, int *select, 
	int *n, double *t, int *ldt, double *vl, int *
	ldvl, double *vr, int *ldvr, int *mm, int *m, 
	double *work, int *info);
 
int dtrexc_(char *compq, int *n, double *t, int *
	ldt, double *q, int *ldq, int *ifst, int *ilst, 
	double *work, int *info);
 
int dtrrfs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, double *a, int *lda, double *b, int *
	ldb, double *x, int *ldx, double *ferr, double *berr, 
	double *work, int *iwork, int *info);
 
int dtrsen_(char *job, char *compq, int *select, int 
	*n, double *t, int *ldt, double *q, int *ldq, 
	double *wr, double *wi, int *m, double *s, double 
	*sep, double *work, int *lwork, int *iwork, int *
	liwork, int *info);
 
int dtrsna_(char *job, char *howmny, int *select, 
	int *n, double *t, int *ldt, double *vl, int *
	ldvl, double *vr, int *ldvr, double *s, double *sep, 
	int *mm, int *m, double *work, int *ldwork, int *
	iwork, int *info);
 
int dtrsyl_(char *trana, char *tranb, int *isgn, int 
	*m, int *n, double *a, int *lda, double *b, int *
	ldb, double *c__, int *ldc, double *scale, int *info);
 
int dtrti2_(char *uplo, char *diag, int *n, double *
	a, int *lda, int *info);
 
int dtrtri_(char *uplo, char *diag, int *n, double *
	a, int *lda, int *info);
 
int dtrtrs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, double *a, int *lda, double *b, int *
	ldb, int *info);
 
int dtzrqf_(int *m, int *n, double *a, int *
	lda, double *tau, int *info);
 
int dtzrzf_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *lwork, int *info);
 
int icmax1_(int *n, single_cx *cx, int *incx);
 
int ieeeck_(int *ispec, float *zero, float *one);
 
int ilaenv_(int *ispec, char *name__, char *opts, int *n1, 
	int *n2, int *n3, int *n4, short name_len, short 
	opts_len);
 
int izmax1_(int *n, double_cx *cx, int *incx);
 
int sbdsdc_(char *uplo, char *compq, int *n, float *d__, 
	float *e, float *u, int *ldu, float *vt, int *ldvt, float *q, 
	int *iq, float *work, int *iwork, int *info);
 
int sbdsqr_(char *uplo, int *n, int *ncvt, int *
	nru, int *ncc, float *d__, float *e, float *vt, int *ldvt, float *
	u, int *ldu, float *c__, int *ldc, float *work, int *info);
 
int sdisna_(char *job, int *m, int *n, float *d__, 
	float *sep, int *info);
 
int sgbbrd_(char *vect, int *m, int *n, int *ncc,
	 int *kl, int *ku, float *ab, int *ldab, float *d__, float *
	e, float *q, int *ldq, float *pt, int *ldpt, float *c__, int 
	*ldc, float *work, int *info);
 
int sgbcon_(char *norm, int *n, int *kl, int *ku,
	 float *ab, int *ldab, int *ipiv, float *anorm, float *rcond, 
	float *work, int *iwork, int *info);
 
int sgbequ_(int *m, int *n, int *kl, int *ku,
	 float *ab, int *ldab, float *r__, float *c__, float *rowcnd, float *
	colcnd, float *amax, int *info);
 
int sgbrfs_(char *trans, int *n, int *kl, int *
	ku, int *nrhs, float *ab, int *ldab, float *afb, int *ldafb,
	 int *ipiv, float *b, int *ldb, float *x, int *ldx, float *
	ferr, float *berr, float *work, int *iwork, int *info);
 
int sgbsv_(int *n, int *kl, int *ku, int *
	nrhs, float *ab, int *ldab, int *ipiv, float *b, int *ldb, 
	int *info);
 
int sgbsvx_(char *fact, char *trans, int *n, int *kl,
	 int *ku, int *nrhs, float *ab, int *ldab, float *afb, 
	int *ldafb, int *ipiv, char *equed, float *r__, float *c__, 
	float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr,
	 float *berr, float *work, int *iwork, int *info);
 
int sgbtf2_(int *m, int *n, int *kl, int *ku,
	 float *ab, int *ldab, int *ipiv, int *info);
 
int sgbtrf_(int *m, int *n, int *kl, int *ku,
	 float *ab, int *ldab, int *ipiv, int *info);
 
int sgbtrs_(char *trans, int *n, int *kl, int *
	ku, int *nrhs, float *ab, int *ldab, int *ipiv, float *b, 
	int *ldb, int *info);
 
int sgebak_(char *job, char *side, int *n, int *ilo, 
	int *ihi, float *scale, int *m, float *v, int *ldv, int 
	*info);
 
int sgebal_(char *job, int *n, float *a, int *lda, 
	int *ilo, int *ihi, float *scale, int *info);
 
int sgebd2_(int *m, int *n, float *a, int *lda, 
	float *d__, float *e, float *tauq, float *taup, float *work, int *info);
 
int sgebrd_(int *m, int *n, float *a, int *lda, 
	float *d__, float *e, float *tauq, float *taup, float *work, int *
	lwork, int *info);
 
int sgecon_(char *norm, int *n, float *a, int *lda, 
	float *anorm, float *rcond, float *work, int *iwork, int *info);
 
int sgeequ_(int *m, int *n, float *a, int *lda, 
	float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, int 
	*info);
 
int sgees_(char *jobvs, char *sort, L_fp select, int *n, 
	float *a, int *lda, int *sdim, float *wr, float *wi, float *vs, 
	int *ldvs, float *work, int *lwork, int *bwork, int *
	info);
 
int sgeesx_(char *jobvs, char *sort, L_fp select, char *
	sense, int *n, float *a, int *lda, int *sdim, float *wr, 
	float *wi, float *vs, int *ldvs, float *rconde, float *rcondv, float *
	work, int *lwork, int *iwork, int *liwork, int *bwork,
	 int *info);
 
int sgeev_(char *jobvl, char *jobvr, int *n, float *a, 
	int *lda, float *wr, float *wi, float *vl, int *ldvl, float *vr, 
	int *ldvr, float *work, int *lwork, int *info);
 
int sgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, int *n, float *a, int *lda, float *wr, float *wi, float *
	vl, int *ldvl, float *vr, int *ldvr, int *ilo, int *
	ihi, float *scale, float *abnrm, float *rconde, float *rcondv, float *work,
	 int *lwork, int *iwork, int *info);
 
int sgegs_(char *jobvsl, char *jobvsr, int *n, float *a, 
	int *lda, float *b, int *ldb, float *alphar, float *alphai, float 
	*beta, float *vsl, int *ldvsl, float *vsr, int *ldvsr, float *
	work, int *lwork, int *info);
 
int sgegv_(char *jobvl, char *jobvr, int *n, float *a, 
	int *lda, float *b, int *ldb, float *alphar, float *alphai, float 
	*beta, float *vl, int *ldvl, float *vr, int *ldvr, float *work, 
	int *lwork, int *info);
 
int sgehd2_(int *n, int *ilo, int *ihi, float *a, 
	int *lda, float *tau, float *work, int *info);
 
int sgehrd_(int *n, int *ilo, int *ihi, float *a, 
	int *lda, float *tau, float *work, int *lwork, int *info);
 
int sgelq2_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *info);
 
int sgelqf_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info);
 
int sgels_(char *trans, int *m, int *n, int *
	nrhs, float *a, int *lda, float *b, int *ldb, float *work, 
	int *lwork, int *info);
 
int sgelsd_(int *m, int *n, int *nrhs, float *a, 
	int *lda, float *b, int *ldb, float *s, float *rcond, int *
	rank, float *work, int *lwork, int *iwork, int *info);
 
int sgelss_(int *m, int *n, int *nrhs, float *a, 
	int *lda, float *b, int *ldb, float *s, float *rcond, int *
	rank, float *work, int *lwork, int *info);
 
int sgelsx_(int *m, int *n, int *nrhs, float *a, 
	int *lda, float *b, int *ldb, int *jpvt, float *rcond, 
	int *rank, float *work, int *info);
 
int sgelsy_(int *m, int *n, int *nrhs, float *a, 
	int *lda, float *b, int *ldb, int *jpvt, float *rcond, 
	int *rank, float *work, int *lwork, int *info);
 
int sgeql2_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *info);
 
int sgeqlf_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info);
 
int sgeqp3_(int *m, int *n, float *a, int *lda, 
	int *jpvt, float *tau, float *work, int *lwork, int *info);
 
int sgeqpf_(int *m, int *n, float *a, int *lda, 
	int *jpvt, float *tau, float *work, int *info);
 
int sgeqr2_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *info);
 
int sgeqrf_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info);
 
int sgerfs_(char *trans, int *n, int *nrhs, float *a, 
	int *lda, float *af, int *ldaf, int *ipiv, float *b, 
	int *ldb, float *x, int *ldx, float *ferr, float *berr, float *
	work, int *iwork, int *info);
 
int sgerq2_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *info);
 
int sgerqf_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info);
 
int sgesc2_(int *n, float *a, int *lda, float *rhs, 
	int *ipiv, int *jpiv, float *scale);
 
int sgesdd_(char *jobz, int *m, int *n, float *a, 
	int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt,
	 float *work, int *lwork, int *iwork, int *info);
 
int sgesv_(int *n, int *nrhs, float *a, int *lda, 
	int *ipiv, float *b, int *ldb, int *info);
 
int sgesvd_(char *jobu, char *jobvt, int *m, int *n, 
	float *a, int *lda, float *s, float *u, int *ldu, float *vt, 
	int *ldvt, float *work, int *lwork, int *info);
 
int sgesvx_(char *fact, char *trans, int *n, int *
	nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, 
	char *equed, float *r__, float *c__, float *b, int *ldb, float *x, 
	int *ldx, float *rcond, float *ferr, float *berr, float *work, 
	int *iwork, int *info);
 
int sgetc2_(int *n, float *a, int *lda, int *ipiv,
	 int *jpiv, int *info);
 
int sgetf2_(int *m, int *n, float *a, int *lda, 
	int *ipiv, int *info);
 
int sgetrf_(int *m, int *n, float *a, int *lda, 
	int *ipiv, int *info);
 
int sgetri_(int *n, float *a, int *lda, int *ipiv,
	 float *work, int *lwork, int *info);
 
int sgetrs_(char *trans, int *n, int *nrhs, float *a, 
	int *lda, int *ipiv, float *b, int *ldb, int *info);
 
int sggbak_(char *job, char *side, int *n, int *ilo, 
	int *ihi, float *lscale, float *rscale, int *m, float *v, 
	int *ldv, int *info);
 
int sggbal_(char *job, int *n, float *a, int *lda, 
	float *b, int *ldb, int *ilo, int *ihi, float *lscale, float 
	*rscale, float *work, int *info);
 
int sgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, int *n, float *a, int *lda, float *b, int *ldb, 
	int *sdim, float *alphar, float *alphai, float *beta, float *vsl, 
	int *ldvsl, float *vsr, int *ldvsr, float *work, int *lwork,
	 int *bwork, int *info);
 
int sggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, char *sense, int *n, float *a, int *lda, float *b, 
	int *ldb, int *sdim, float *alphar, float *alphai, float *beta, 
	float *vsl, int *ldvsl, float *vsr, int *ldvsr, float *rconde, 
	float *rcondv, float *work, int *lwork, int *iwork, int *
	liwork, int *bwork, int *info);
 
int sggev_(char *jobvl, char *jobvr, int *n, float *a, 
	int *lda, float *b, int *ldb, float *alphar, float *alphai, float 
	*beta, float *vl, int *ldvl, float *vr, int *ldvr, float *work, 
	int *lwork, int *info);
 
int sggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, int *n, float *a, int *lda, float *b, int *ldb, float 
	*alphar, float *alphai, float *beta, float *vl, int *ldvl, float *vr, 
	int *ldvr, int *ilo, int *ihi, float *lscale, float *rscale,
	 float *abnrm, float *bbnrm, float *rconde, float *rcondv, float *work, 
	int *lwork, int *iwork, int *bwork, int *info);
 
int sggglm_(int *n, int *m, int *p, float *a, 
	int *lda, float *b, int *ldb, float *d__, float *x, float *y, 
	float *work, int *lwork, int *info);
 
int sgghrd_(char *compq, char *compz, int *n, int *
	ilo, int *ihi, float *a, int *lda, float *b, int *ldb, float 
	*q, int *ldq, float *z__, int *ldz, int *info);
 
int sgglse_(int *m, int *n, int *p, float *a, 
	int *lda, float *b, int *ldb, float *c__, float *d__, float *x, 
	float *work, int *lwork, int *info);
 
int sggqrf_(int *n, int *m, int *p, float *a, 
	int *lda, float *taua, float *b, int *ldb, float *taub, float *
	work, int *lwork, int *info);
 
int sggrqf_(int *m, int *p, int *n, float *a, 
	int *lda, float *taua, float *b, int *ldb, float *taub, float *
	work, int *lwork, int *info);
 
int sggsvd_(char *jobu, char *jobv, char *jobq, int *m, 
	int *n, int *p, int *k, int *l, float *a, int *lda,
	 float *b, int *ldb, float *alpha, float *beta, float *u, int *
	ldu, float *v, int *ldv, float *q, int *ldq, float *work, 
	int *iwork, int *info);
 
int sggsvp_(char *jobu, char *jobv, char *jobq, int *m, 
	int *p, int *n, float *a, int *lda, float *b, int *ldb, 
	float *tola, float *tolb, int *k, int *l, float *u, int *ldu,
	 float *v, int *ldv, float *q, int *ldq, int *iwork, float *
	tau, float *work, int *info);
 
int sgtcon_(char *norm, int *n, float *dl, float *d__, 
	float *du, float *du2, int *ipiv, float *anorm, float *rcond, float *
	work, int *iwork, int *info);
 
int sgtrfs_(char *trans, int *n, int *nrhs, float *dl,
	 float *d__, float *du, float *dlf, float *df, float *duf, float *du2, 
	int *ipiv, float *b, int *ldb, float *x, int *ldx, float *
	ferr, float *berr, float *work, int *iwork, int *info);
 
int sgtsv_(int *n, int *nrhs, float *dl, float *d__, 
	float *du, float *b, int *ldb, int *info);
 
int sgtsvx_(char *fact, char *trans, int *n, int *
	nrhs, float *dl, float *d__, float *du, float *dlf, float *df, float *duf, 
	float *du2, int *ipiv, float *b, int *ldb, float *x, int *
	ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, 
	int *info);
 
int sgttrf_(int *n, float *dl, float *d__, float *du, float *
	du2, int *ipiv, int *info);
 
int sgttrs_(char *trans, int *n, int *nrhs, float *dl,
	 float *d__, float *du, float *du2, int *ipiv, float *b, int *ldb,
	 int *info);
 
int sgtts2_(int *itrans, int *n, int *nrhs, float 
	*dl, float *d__, float *du, float *du2, int *ipiv, float *b, int *
	ldb);
 
int shgeqz_(char *job, char *compq, char *compz, int *n, 
	int *ilo, int *ihi, float *a, int *lda, float *b, int *
	ldb, float *alphar, float *alphai, float *beta, float *q, int *ldq, 
	float *z__, int *ldz, float *work, int *lwork, int *info);
 
int shsein_(char *side, char *eigsrc, char *initv, int *
	select, int *n, float *h__, int *ldh, float *wr, float *wi, float 
	*vl, int *ldvl, float *vr, int *ldvr, int *mm, int *m, 
	float *work, int *ifaill, int *ifailr, int *info);
 
int shseqr_(char *job, char *compz, int *n, int *ilo,
	 int *ihi, float *h__, int *ldh, float *wr, float *wi, float *z__,
	 int *ldz, float *work, int *lwork, int *info);
 
int slabad_(float *small, float *large);
 
int slabrd_(int *m, int *n, int *nb, float *a, 
	int *lda, float *d__, float *e, float *tauq, float *taup, float *x, 
	int *ldx, float *y, int *ldy);
 
int slacon_(int *n, float *v, float *x, int *isgn, 
	float *est, int *kase);
 
int slacpy_(char *uplo, int *m, int *n, float *a, 
	int *lda, float *b, int *ldb);
 
int sladiv_(float *a, float *b, float *c__, float *d__, float *p, 
	float *q);
 
int slae2_(float *a, float *b, float *c__, float *rt1, float *rt2);
 
int slaebz_(int *ijob, int *nitmax, int *n, 
	int *mmax, int *minp, int *nbmin, float *abstol, float *
	reltol, float *pivmin, float *d__, float *e, float *e2, int *nval, 
	float *ab, float *c__, int *mout, int *nab, float *work, int 
	*iwork, int *info);
 
int slaed0_(int *icompq, int *qsiz, int *n, float 
	*d__, float *e, float *q, int *ldq, float *qstore, int *ldqs, 
	float *work, int *iwork, int *info);
 
int slaed1_(int *n, float *d__, float *q, int *ldq, 
	int *indxq, float *rho, int *cutpnt, float *work, int *
	iwork, int *info);
 
int slaed2_(int *k, int *n, int *n1, float *d__, 
	float *q, int *ldq, int *indxq, float *rho, float *z__, float *
	dlamda, float *w, float *q2, int *indx, int *indxc, int *
	indxp, int *coltyp, int *info);
 
int slaed3_(int *k, int *n, int *n1, float *d__, 
	float *q, int *ldq, float *rho, float *dlamda, float *q2, int *
	indx, int *ctot, float *w, float *s, int *info);
 
int slaed4_(int *n, int *i__, float *d__, float *z__, 
	float *delta, float *rho, float *dlam, int *info);
 
int slaed5_(int *i__, float *d__, float *z__, float *delta, 
	float *rho, float *dlam);
 
int slaed6_(int *kniter, int *orgati, float *rho, 
	float *d__, float *z__, float *finit, float *tau, int *info);
 
int slaed7_(int *icompq, int *n, int *qsiz, 
	int *tlvls, int *curlvl, int *curpbm, float *d__, float *q, 
	int *ldq, int *indxq, float *rho, int *cutpnt, float *
	qstore, int *qptr, int *prmptr, int *perm, int *
	givptr, int *givcol, float *givnum, float *work, int *iwork, 
	int *info);
 
int slaed8_(int *icompq, int *k, int *n, int 
	*qsiz, float *d__, float *q, int *ldq, int *indxq, float *rho, 
	int *cutpnt, float *z__, float *dlamda, float *q2, int *ldq2, 
	float *w, int *perm, int *givptr, int *givcol, float *
	givnum, int *indxp, int *indx, int *info);
 
int slaed9_(int *k, int *kstart, int *kstop, 
	int *n, float *d__, float *q, int *ldq, float *rho, float *dlamda,
	 float *w, float *s, int *lds, int *info);
 
int slaeda_(int *n, int *tlvls, int *curlvl, 
	int *curpbm, int *prmptr, int *perm, int *givptr, 
	int *givcol, float *givnum, float *q, int *qptr, float *z__, 
	float *ztemp, int *info);
 
int slaein_(int *rightv, int *noinit, int *n, 
	float *h__, int *ldh, float *wr, float *wi, float *vr, float *vi, float 
	*b, int *ldb, float *work, float *eps3, float *smlnum, float *bignum, 
	int *info);
 
int slaev2_(float *a, float *b, float *c__, float *rt1, float *
	rt2, float *cs1, float *sn1);
 
int slaexc_(int *wantq, int *n, float *t, int *
	ldt, float *q, int *ldq, int *j1, int *n1, int *n2, 
	float *work, int *info);
 
int slag2_(float *a, int *lda, float *b, int *ldb, 
	float *safmin, float *scale1, float *scale2, float *wr1, float *wr2, float *
	wi);
 
int slags2_(int *upper, float *a1, float *a2, float *a3, 
	float *b1, float *b2, float *b3, float *csu, float *snu, float *csv, float *
	snv, float *csq, float *snq);
 
int slagtf_(int *n, float *a, float *lambda, float *b, float 
	*c__, float *tol, float *d__, int *in, int *info);
 
int slagtm_(char *trans, int *n, int *nrhs, float *
	alpha, float *dl, float *d__, float *du, float *x, int *ldx, float *
	beta, float *b, int *ldb);
 
int slagts_(int *job, int *n, float *a, float *b, float 
	*c__, float *d__, int *in, float *y, float *tol, int *info);
 
int slagv2_(float *a, int *lda, float *b, int *ldb, 
	float *alphar, float *alphai, float *beta, float *csl, float *snl, float *
	csr, float *snr);
 
int slahqr_(int *wantt, int *wantz, int *n, 
	int *ilo, int *ihi, float *h__, int *ldh, float *wr, float *
	wi, int *iloz, int *ihiz, float *z__, int *ldz, int *
	info);
 
int slahrd_(int *n, int *k, int *nb, float *a, 
	int *lda, float *tau, float *t, int *ldt, float *y, int *ldy);
 
int slaic1_(int *job, int *j, float *x, float *sest, 
	float *w, float *gamma, float *sestpr, float *s, float *c__);
 
int slaln2_(int *ltrans, int *na, int *nw, float *
	smin, float *ca, float *a, int *lda, float *d1, float *d2, float *b, 
	int *ldb, float *wr, float *wi, float *x, int *ldx, float *scale, 
	float *xnorm, int *info);
 
int slals0_(int *icompq, int *nl, int *nr, 
	int *sqre, int *nrhs, float *b, int *ldb, float *bx, 
	int *ldbx, int *perm, int *givptr, int *givcol, 
	int *ldgcol, float *givnum, int *ldgnum, float *poles, float *
	difl, float *difr, float *z__, int *k, float *c__, float *s, float *
	work, int *info);
 
int slalsa_(int *icompq, int *smlsiz, int *n, 
	int *nrhs, float *b, int *ldb, float *bx, int *ldbx, float *
	u, int *ldu, float *vt, int *k, float *difl, float *difr, float *
	z__, float *poles, int *givptr, int *givcol, int *ldgcol, 
	int *perm, float *givnum, float *c__, float *s, float *work, int *
	iwork, int *info);
 
int slalsd_(char *uplo, int *smlsiz, int *n, int 
	*nrhs, float *d__, float *e, float *b, int *ldb, float *rcond, 
	int *rank, float *work, int *iwork, int *info);
 
int slamc1_(int *beta, int *t, int *rnd, int 
	*ieee1);
 
int slamc2_(int *beta, int *t, int *rnd, float *
	eps, int *emin, float *rmin, int *emax, float *rmax);
 
int slamc4_(int *emin, float *start, int *base);
 
int slamc5_(int *beta, int *p, int *emin, 
	int *ieee, int *emax, float *rmax);
 
int slamrg_(int *n1, int *n2, float *a, int *
	strd1, int *strd2, int *index);
 
int slanv2_(float *a, float *b, float *c__, float *d__, float *
	rt1r, float *rt1i, float *rt2r, float *rt2i, float *cs, float *sn);
 
int slapll_(int *n, float *x, int *incx, float *y, 
	int *incy, float *ssmin);
 
int slapmt_(int *forwrd, int *m, int *n, float *x,
	 int *ldx, int *k);
 
int slaqgb_(int *m, int *n, int *kl, int *ku,
	 float *ab, int *ldab, float *r__, float *c__, float *rowcnd, float *
	colcnd, float *amax, char *equed);
 
int slaqge_(int *m, int *n, float *a, int *lda, 
	float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, char *
	equed);
 
int slaqp2_(int *m, int *n, int *offset, float *a,
	 int *lda, int *jpvt, float *tau, float *vn1, float *vn2, float *
	work);
 
int slaqps_(int *m, int *n, int *offset, int 
	*nb, int *kb, float *a, int *lda, int *jpvt, float *tau, 
	float *vn1, float *vn2, float *auxv, float *f, int *ldf);
 
int slaqsb_(char *uplo, int *n, int *kd, float *ab, 
	int *ldab, float *s, float *scond, float *amax, char *equed);
 
int slaqsp_(char *uplo, int *n, float *ap, float *s, float *
	scond, float *amax, char *equed);
 
int slaqsy_(char *uplo, int *n, float *a, int *lda, 
	float *s, float *scond, float *amax, char *equed);
 
int slaqtr_(int *ltran, int *lfloat, int *n, float 
	*t, int *ldt, float *b, float *w, float *scale, float *x, float *work, 
	int *info);
 
int slar1v_(int *n, int *b1, int *bn, float *
	sigma, float *d__, float *l, float *ld, float *lld, float *gersch, float *
	z__, float *ztz, float *mingma, int *r__, int *isuppz, float *
	work);
 
int slar2v_(int *n, float *x, float *y, float *z__, int 
	*incx, float *c__, float *s, int *incc);
 
int slarf_(char *side, int *m, int *n, float *v, 
	int *incv, float *tau, float *c__, int *ldc, float *work);
 
int slarfb_(char *side, char *trans, char *direct, char *
	storev, int *m, int *n, int *k, float *v, int *ldv, 
	float *t, int *ldt, float *c__, int *ldc, float *work, int *
	ldwork);
 
int slarfg_(int *n, float *alpha, float *x, int *incx, 
	float *tau);
 
int slarft_(char *direct, char *storev, int *n, int *
	k, float *v, int *ldv, float *tau, float *t, int *ldt);
 
int slarfx_(char *side, int *m, int *n, float *v, 
	float *tau, float *c__, int *ldc, float *work);
 
int slargv_(int *n, float *x, int *incx, float *y, 
	int *incy, float *c__, int *incc);
 
int slarnv_(int *idist, int *iseed, int *n, float 
	*x);
 
int slarrb_(int *n, float *d__, float *l, float *ld, float *
	lld, int *ifirst, int *ilast, float *sigma, float *reltol, float 
	*w, float *wgap, float *werr, float *work, int *iwork, int *info);
 
int slarre_(int *n, float *d__, float *e, float *tol, 
	int *nsplit, int *isplit, int *m, float *w, float *woff, 
	float *gersch, float *work, int *info);
 
int slarrf_(int *n, float *d__, float *l, float *ld, float *
	lld, int *ifirst, int *ilast, float *w, float *dplus, float *
	lplus, float *work, int *iwork, int *info);
 
int slarrv_(int *n, float *d__, float *l, int *isplit, 
	int *m, float *w, int *iblock, float *gersch, float *tol, float *
	z__, int *ldz, int *isuppz, float *work, int *iwork, 
	int *info);
 
int slartg_(float *f, float *g, float *cs, float *sn, float *r__);
 
int slartv_(int *n, float *x, int *incx, float *y, 
	int *incy, float *c__, float *s, int *incc);
 
int slaruv_(int *iseed, int *n, float *x);
 
int slarz_(char *side, int *m, int *n, int *l, 
	float *v, int *incv, float *tau, float *c__, int *ldc, float *
	work);
 
int slarzb_(char *side, char *trans, char *direct, char *
	storev, int *m, int *n, int *k, int *l, float *v, 
	int *ldv, float *t, int *ldt, float *c__, int *ldc, float *
	work, int *ldwork);
 
int slarzt_(char *direct, char *storev, int *n, int *
	k, float *v, int *ldv, float *tau, float *t, int *ldt);
 
int slas2_(float *f, float *g, float *h__, float *ssmin, float *
	ssmax);
 
int slascl_(char *type__, int *kl, int *ku, float *
	cfrom, float *cto, int *m, int *n, float *a, int *lda, 
	int *info);
 
int slasd0_(int *n, int *sqre, float *d__, float *e, 
	float *u, int *ldu, float *vt, int *ldvt, int *smlsiz, 
	int *iwork, float *work, int *info);
 
int slasd1_(int *nl, int *nr, int *sqre, float *
	d__, float *alpha, float *beta, float *u, int *ldu, float *vt, 
	int *ldvt, int *idxq, int *iwork, float *work, int *
	info);
 
int slasd2_(int *nl, int *nr, int *sqre, int 
	*k, float *d__, float *z__, float *alpha, float *beta, float *u, int *
	ldu, float *vt, int *ldvt, float *dsigma, float *u2, int *ldu2, 
	float *vt2, int *ldvt2, int *idxp, int *idx, int *idxc,
	 int *idxq, int *coltyp, int *info);
 
int slasd3_(int *nl, int *nr, int *sqre, int 
	*k, float *d__, float *q, int *ldq, float *dsigma, float *u, int *
	ldu, float *u2, int *ldu2, float *vt, int *ldvt, float *vt2, 
	int *ldvt2, int *idxc, int *ctot, float *z__, int *
	info);
 
int slasd4_(int *n, int *i__, float *d__, float *z__, 
	float *delta, float *rho, float *sigma, float *work, int *info);
 
int slasd5_(int *i__, float *d__, float *z__, float *delta, 
	float *rho, float *dsigma, float *work);
 
int slasd6_(int *icompq, int *nl, int *nr, 
	int *sqre, float *d__, float *vf, float *vl, float *alpha, float *beta,
	 int *idxq, int *perm, int *givptr, int *givcol, 
	int *ldgcol, float *givnum, int *ldgnum, float *poles, float *
	difl, float *difr, float *z__, int *k, float *c__, float *s, float *
	work, int *iwork, int *info);
 
int slasd7_(int *icompq, int *nl, int *nr, 
	int *sqre, int *k, float *d__, float *z__, float *zw, float *vf, 
	float *vfw, float *vl, float *vlw, float *alpha, float *beta, float *dsigma,
	 int *idx, int *idxp, int *idxq, int *perm, int *
	givptr, int *givcol, int *ldgcol, float *givnum, int *
	ldgnum, float *c__, float *s, int *info);
 
int slasd8_(int *icompq, int *k, float *d__, float *
	z__, float *vf, float *vl, float *difl, float *difr, int *lddifr, 
	float *dsigma, float *work, int *info);
 
int slasd9_(int *icompq, int *ldu, int *k, float *
	d__, float *z__, float *vf, float *vl, float *difl, float *difr, float *
	dsigma, float *work, int *info);
 
int slasda_(int *icompq, int *smlsiz, int *n, 
	int *sqre, float *d__, float *e, float *u, int *ldu, float *vt, 
	int *k, float *difl, float *difr, float *z__, float *poles, int *
	givptr, int *givcol, int *ldgcol, int *perm, float *givnum,
	 float *c__, float *s, float *work, int *iwork, int *info);
 
int slasdq_(char *uplo, int *sqre, int *n, int *
	ncvt, int *nru, int *ncc, float *d__, float *e, float *vt, 
	int *ldvt, float *u, int *ldu, float *c__, int *ldc, float *
	work, int *info);
 
int slasdt_(int *n, int *lvl, int *nd, int *
	inode, int *ndiml, int *ndimr, int *msub);
 
int slaset_(char *uplo, int *m, int *n, float *alpha, 
	float *beta, float *a, int *lda);
 
int slasq1_(int *n, float *d__, float *e, float *work, 
	int *info);
 
int slasq2_(int *n, float *z__, int *info);
 
int slasq3_(int *i0, int *n0, float *z__, int *pp,
	 float *dmin__, float *sigma, float *desig, float *qmax, int *nfail, 
	int *iter, int *ndiv, int *ieee);
 
int slasq4_(int *i0, int *n0, float *z__, int *pp,
	 int *n0in, float *dmin__, float *dmin1, float *dmin2, float *dn, 
	float *dn1, float *dn2, float *tau, int *ttype);
 
int slasq5_(int *i0, int *n0, float *z__, int *pp,
	 float *tau, float *dmin__, float *dmin1, float *dmin2, float *dn, float *
	dnm1, float *dnm2, int *ieee);
 
int slasq6_(int *i0, int *n0, float *z__, int *pp,
	 float *dmin__, float *dmin1, float *dmin2, float *dn, float *dnm1, float *
	dnm2);
 
int slasr_(char *side, char *pivot, char *direct, int *m,
	 int *n, float *c__, float *s, float *a, int *lda);
 
int slasrt_(char *id, int *n, float *d__, int *info);
 
int slassq_(int *n, float *x, int *incx, float *scale, 
	float *sumsq);
 
int slasv2_(float *f, float *g, float *h__, float *ssmin, float *
	ssmax, float *snr, float *csr, float *snl, float *csl);
 
int slaswp_(int *n, float *a, int *lda, int *k1, 
	int *k2, int *ipiv, int *incx);
 
int slasy2_(int *ltranl, int *ltranr, int *isgn, 
	int *n1, int *n2, float *tl, int *ldtl, float *tr, int *
	ldtr, float *b, int *ldb, float *scale, float *x, int *ldx, float 
	*xnorm, int *info);
 
int slasyf_(char *uplo, int *n, int *nb, int *kb,
	 float *a, int *lda, int *ipiv, float *w, int *ldw, int 
	*info);
 
int slatbs_(char *uplo, char *trans, char *diag, char *
	normin, int *n, int *kd, float *ab, int *ldab, float *x, 
	float *scale, float *cnorm, int *info);
 
int slatdf_(int *ijob, int *n, float *z__, int *
	ldz, float *rhs, float *rdsum, float *rdscal, int *ipiv, int *
	jpiv);
 
int slatps_(char *uplo, char *trans, char *diag, char *
	normin, int *n, float *ap, float *x, float *scale, float *cnorm, 
	int *info);
 
int slatrd_(char *uplo, int *n, int *nb, float *a, 
	int *lda, float *e, float *tau, float *w, int *ldw);
 
int slatrs_(char *uplo, char *trans, char *diag, char *
	normin, int *n, float *a, int *lda, float *x, float *scale, float 
	*cnorm, int *info);
 
int slatrz_(int *m, int *n, int *l, float *a, 
	int *lda, float *tau, float *work);
 
int slatzm_(char *side, int *m, int *n, float *v, 
	int *incv, float *tau, float *c1, float *c2, int *ldc, float *
	work);
 
int slauu2_(char *uplo, int *n, float *a, int *lda, 
	int *info);
 
int slauum_(char *uplo, int *n, float *a, int *lda, 
	int *info);
 
int sopgtr_(char *uplo, int *n, float *ap, float *tau, 
	float *q, int *ldq, float *work, int *info);
 
int sopmtr_(char *side, char *uplo, char *trans, int *m, 
	int *n, float *ap, float *tau, float *c__, int *ldc, float *work, 
	int *info);
 
int sorg2l_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *info);
 
int sorg2r_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *info);
 
int sorgbr_(char *vect, int *m, int *n, int *k, 
	float *a, int *lda, float *tau, float *work, int *lwork, int 
	*info);
 
int sorghr_(int *n, int *ilo, int *ihi, float *a, 
	int *lda, float *tau, float *work, int *lwork, int *info);
 
int sorgl2_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *info);
 
int sorglq_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *lwork, int *info);
 
int sorgql_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *lwork, int *info);
 
int sorgqr_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *lwork, int *info);
 
int sorgr2_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *info);
 
int sorgrq_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *lwork, int *info);
 
int sorgtr_(char *uplo, int *n, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info);
 
int sorm2l_(char *side, char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *info);
 
int sorm2r_(char *side, char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *info);
 
int sormbr_(char *vect, char *side, char *trans, int *m, 
	int *n, int *k, float *a, int *lda, float *tau, float *c__, 
	int *ldc, float *work, int *lwork, int *info);
 
int sormhr_(char *side, char *trans, int *m, int *n, 
	int *ilo, int *ihi, float *a, int *lda, float *tau, float *
	c__, int *ldc, float *work, int *lwork, int *info);
 
int sorml2_(char *side, char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *info);
 
int sormlq_(char *side, char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *lwork, int *info);
 
int sormql_(char *side, char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *lwork, int *info);
 
int sormqr_(char *side, char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *lwork, int *info);
 
int sormr2_(char *side, char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *info);
 
int sormr3_(char *side, char *trans, int *m, int *n, 
	int *k, int *l, float *a, int *lda, float *tau, float *c__, 
	int *ldc, float *work, int *info);
 
int sormrq_(char *side, char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *lwork, int *info);
 
int sormrz_(char *side, char *trans, int *m, int *n, 
	int *k, int *l, float *a, int *lda, float *tau, float *c__, 
	int *ldc, float *work, int *lwork, int *info);
 
int sormtr_(char *side, char *uplo, char *trans, int *m, 
	int *n, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *lwork, int *info);
 
int spbcon_(char *uplo, int *n, int *kd, float *ab, 
	int *ldab, float *anorm, float *rcond, float *work, int *iwork, 
	int *info);
 
int spbequ_(char *uplo, int *n, int *kd, float *ab, 
	int *ldab, float *s, float *scond, float *amax, int *info);
 
int spbrfs_(char *uplo, int *n, int *kd, int *
	nrhs, float *ab, int *ldab, float *afb, int *ldafb, float *b, 
	int *ldb, float *x, int *ldx, float *ferr, float *berr, float *
	work, int *iwork, int *info);
 
int spbstf_(char *uplo, int *n, int *kd, float *ab, 
	int *ldab, int *info);
 
int spbsv_(char *uplo, int *n, int *kd, int *
	nrhs, float *ab, int *ldab, float *b, int *ldb, int *info);
 
int spbsvx_(char *fact, char *uplo, int *n, int *kd, 
	int *nrhs, float *ab, int *ldab, float *afb, int *ldafb, 
	char *equed, float *s, float *b, int *ldb, float *x, int *ldx, 
	float *rcond, float *ferr, float *berr, float *work, int *iwork, 
	int *info);
 
int spbtf2_(char *uplo, int *n, int *kd, float *ab, 
	int *ldab, int *info);
 
int spbtrf_(char *uplo, int *n, int *kd, float *ab, 
	int *ldab, int *info);
 
int spbtrs_(char *uplo, int *n, int *kd, int *
	nrhs, float *ab, int *ldab, float *b, int *ldb, int *info);
 
int spocon_(char *uplo, int *n, float *a, int *lda, 
	float *anorm, float *rcond, float *work, int *iwork, int *info);
 
int spoequ_(int *n, float *a, int *lda, float *s, float 
	*scond, float *amax, int *info);
 
int sporfs_(char *uplo, int *n, int *nrhs, float *a, 
	int *lda, float *af, int *ldaf, float *b, int *ldb, float *x,
	 int *ldx, float *ferr, float *berr, float *work, int *iwork, 
	int *info);
 
int sposv_(char *uplo, int *n, int *nrhs, float *a, 
	int *lda, float *b, int *ldb, int *info);
 
int sposvx_(char *fact, char *uplo, int *n, int *
	nrhs, float *a, int *lda, float *af, int *ldaf, char *equed, 
	float *s, float *b, int *ldb, float *x, int *ldx, float *rcond, 
	float *ferr, float *berr, float *work, int *iwork, int *info);
 
int spotf2_(char *uplo, int *n, float *a, int *lda, 
	int *info);
 
int spotrf_(char *uplo, int *n, float *a, int *lda, 
	int *info);
 
int spotri_(char *uplo, int *n, float *a, int *lda, 
	int *info);
 
int spotrs_(char *uplo, int *n, int *nrhs, float *a, 
	int *lda, float *b, int *ldb, int *info);
 
int sppcon_(char *uplo, int *n, float *ap, float *anorm, 
	float *rcond, float *work, int *iwork, int *info);
 
int sppequ_(char *uplo, int *n, float *ap, float *s, float *
	scond, float *amax, int *info);
 
int spprfs_(char *uplo, int *n, int *nrhs, float *ap, 
	float *afp, float *b, int *ldb, float *x, int *ldx, float *ferr, 
	float *berr, float *work, int *iwork, int *info);
 
int sppsv_(char *uplo, int *n, int *nrhs, float *ap, 
	float *b, int *ldb, int *info);
 
int sppsvx_(char *fact, char *uplo, int *n, int *
	nrhs, float *ap, float *afp, char *equed, float *s, float *b, int *
	ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float 
	*work, int *iwork, int *info);
 
int spptrf_(char *uplo, int *n, float *ap, int *info);
 
int spptri_(char *uplo, int *n, float *ap, int *info);
 
int spptrs_(char *uplo, int *n, int *nrhs, float *ap, 
	float *b, int *ldb, int *info);
 
int sptcon_(int *n, float *d__, float *e, float *anorm, 
	float *rcond, float *work, int *info);
 
int spteqr_(char *compz, int *n, float *d__, float *e, 
	float *z__, int *ldz, float *work, int *info);
 
int sptrfs_(int *n, int *nrhs, float *d__, float *e, 
	float *df, float *ef, float *b, int *ldb, float *x, int *ldx, 
	float *ferr, float *berr, float *work, int *info);
 
int sptsv_(int *n, int *nrhs, float *d__, float *e, 
	float *b, int *ldb, int *info);
 
int sptsvx_(char *fact, int *n, int *nrhs, float *d__,
	 float *e, float *df, float *ef, float *b, int *ldb, float *x, int 
	*ldx, float *rcond, float *ferr, float *berr, float *work, int *info);
 
int spttrf_(int *n, float *d__, float *e, int *info);
 
int spttrs_(int *n, int *nrhs, float *d__, float *e, 
	float *b, int *ldb, int *info);
 
int sptts2_(int *n, int *nrhs, float *d__, float *e, 
	float *b, int *ldb);
 
int srscl_(int *n, float *sa, float *sx, int *incx);
 
int ssbev_(char *jobz, char *uplo, int *n, int *kd, 
	float *ab, int *ldab, float *w, float *z__, int *ldz, float *work,
	 int *info);
 
int ssbevd_(char *jobz, char *uplo, int *n, int *kd, 
	float *ab, int *ldab, float *w, float *z__, int *ldz, float *work,
	 int *lwork, int *iwork, int *liwork, int *info);
 
int ssbevx_(char *jobz, char *range, char *uplo, int *n, 
	int *kd, float *ab, int *ldab, float *q, int *ldq, float *vl,
	 float *vu, int *il, int *iu, float *abstol, int *m, float *
	w, float *z__, int *ldz, float *work, int *iwork, int *
	ifail, int *info);
 
int ssbgst_(char *vect, char *uplo, int *n, int *ka, 
	int *kb, float *ab, int *ldab, float *bb, int *ldbb, float *
	x, int *ldx, float *work, int *info);
 
int ssbgv_(char *jobz, char *uplo, int *n, int *ka, 
	int *kb, float *ab, int *ldab, float *bb, int *ldbb, float *
	w, float *z__, int *ldz, float *work, int *info);
 
int ssbgvd_(char *jobz, char *uplo, int *n, int *ka, 
	int *kb, float *ab, int *ldab, float *bb, int *ldbb, float *
	w, float *z__, int *ldz, float *work, int *lwork, int *
	iwork, int *liwork, int *info);
 
int ssbgvx_(char *jobz, char *range, char *uplo, int *n, 
	int *ka, int *kb, float *ab, int *ldab, float *bb, int *
	ldbb, float *q, int *ldq, float *vl, float *vu, int *il, int 
	*iu, float *abstol, int *m, float *w, float *z__, int *ldz, float 
	*work, int *iwork, int *ifail, int *info);
 
int ssbtrd_(char *vect, char *uplo, int *n, int *kd, 
	float *ab, int *ldab, float *d__, float *e, float *q, int *ldq, 
	float *work, int *info);
 
int sspcon_(char *uplo, int *n, float *ap, int *ipiv, 
	float *anorm, float *rcond, float *work, int *iwork, int *info);
 
int sspev_(char *jobz, char *uplo, int *n, float *ap, 
	float *w, float *z__, int *ldz, float *work, int *info);
 
int sspevd_(char *jobz, char *uplo, int *n, float *ap, 
	float *w, float *z__, int *ldz, float *work, int *lwork, int 
	*iwork, int *liwork, int *info);
 
int sspevx_(char *jobz, char *range, char *uplo, int *n, 
	float *ap, float *vl, float *vu, int *il, int *iu, float *abstol, 
	int *m, float *w, float *z__, int *ldz, float *work, int *
	iwork, int *ifail, int *info);
 
int sspgst_(int *itype, char *uplo, int *n, float *ap,
	 float *bp, int *info);
 
int sspgv_(int *itype, char *jobz, char *uplo, int *
	n, float *ap, float *bp, float *w, float *z__, int *ldz, float *work, 
	int *info);
 
int sspgvd_(int *itype, char *jobz, char *uplo, int *
	n, float *ap, float *bp, float *w, float *z__, int *ldz, float *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
int sspgvx_(int *itype, char *jobz, char *range, char *
	uplo, int *n, float *ap, float *bp, float *vl, float *vu, int *il,
	 int *iu, float *abstol, int *m, float *w, float *z__, int *
	ldz, float *work, int *iwork, int *ifail, int *info);
 
int ssprfs_(char *uplo, int *n, int *nrhs, float *ap, 
	float *afp, int *ipiv, float *b, int *ldb, float *x, int *
	ldx, float *ferr, float *berr, float *work, int *iwork, int *
	info);
 
int sspsv_(char *uplo, int *n, int *nrhs, float *ap, 
	int *ipiv, float *b, int *ldb, int *info);
 
int sspsvx_(char *fact, char *uplo, int *n, int *
	nrhs, float *ap, float *afp, int *ipiv, float *b, int *ldb, float 
	*x, int *ldx, float *rcond, float *ferr, float *berr, float *work, 
	int *iwork, int *info);
 
int ssptrd_(char *uplo, int *n, float *ap, float *d__, 
	float *e, float *tau, int *info);
 
int ssptrf_(char *uplo, int *n, float *ap, int *ipiv, 
	int *info);
 
int ssptri_(char *uplo, int *n, float *ap, int *ipiv, 
	float *work, int *info);
 
int ssptrs_(char *uplo, int *n, int *nrhs, float *ap, 
	int *ipiv, float *b, int *ldb, int *info);
 
int sstebz_(char *range, char *order, int *n, float *vl, 
	float *vu, int *il, int *iu, float *abstol, float *d__, float *e, 
	int *m, int *nsplit, float *w, int *iblock, int *
	isplit, float *work, int *iwork, int *info);
 
int sstedc_(char *compz, int *n, float *d__, float *e, 
	float *z__, int *ldz, float *work, int *lwork, int *iwork, 
	int *liwork, int *info);
 
int sstegr_(char *jobz, char *range, int *n, float *d__, 
	float *e, float *vl, float *vu, int *il, int *iu, float *abstol, 
	int *m, float *w, float *z__, int *ldz, int *isuppz, float *
	work, int *lwork, int *iwork, int *liwork, int *info);
 
int sstein_(int *n, float *d__, float *e, int *m, float 
	*w, int *iblock, int *isplit, float *z__, int *ldz, float *
	work, int *iwork, int *ifail, int *info);
 
int ssteqr_(char *compz, int *n, float *d__, float *e, 
	float *z__, int *ldz, float *work, int *info);
 
int ssterf_(int *n, float *d__, float *e, int *info);
 
int sstev_(char *jobz, int *n, float *d__, float *e, float *
	z__, int *ldz, float *work, int *info);
 
int sstevd_(char *jobz, int *n, float *d__, float *e, float 
	*z__, int *ldz, float *work, int *lwork, int *iwork, 
	int *liwork, int *info);
 
int sstevr_(char *jobz, char *range, int *n, float *d__, 
	float *e, float *vl, float *vu, int *il, int *iu, float *abstol, 
	int *m, float *w, float *z__, int *ldz, int *isuppz, float *
	work, int *lwork, int *iwork, int *liwork, int *info);
 
int sstevx_(char *jobz, char *range, int *n, float *d__, 
	float *e, float *vl, float *vu, int *il, int *iu, float *abstol, 
	int *m, float *w, float *z__, int *ldz, float *work, int *
	iwork, int *ifail, int *info);
 
int ssycon_(char *uplo, int *n, float *a, int *lda, 
	int *ipiv, float *anorm, float *rcond, float *work, int *iwork, 
	int *info);
 
int ssyev_(char *jobz, char *uplo, int *n, float *a, 
	int *lda, float *w, float *work, int *lwork, int *info);
 
int ssyevd_(char *jobz, char *uplo, int *n, float *a, 
	int *lda, float *w, float *work, int *lwork, int *iwork, 
	int *liwork, int *info);
 
int ssyevr_(char *jobz, char *range, char *uplo, int *n, 
	float *a, int *lda, float *vl, float *vu, int *il, int *iu, 
	float *abstol, int *m, float *w, float *z__, int *ldz, int *
	isuppz, float *work, int *lwork, int *iwork, int *liwork, 
	int *info);
 
int ssyevx_(char *jobz, char *range, char *uplo, int *n, 
	float *a, int *lda, float *vl, float *vu, int *il, int *iu, 
	float *abstol, int *m, float *w, float *z__, int *ldz, float *
	work, int *lwork, int *iwork, int *ifail, int *info);
 
int ssygs2_(int *itype, char *uplo, int *n, float *a, 
	int *lda, float *b, int *ldb, int *info);
 
int ssygst_(int *itype, char *uplo, int *n, float *a, 
	int *lda, float *b, int *ldb, int *info);
 
int ssygv_(int *itype, char *jobz, char *uplo, int *
	n, float *a, int *lda, float *b, int *ldb, float *w, float *work, 
	int *lwork, int *info);
 
int ssygvd_(int *itype, char *jobz, char *uplo, int *
	n, float *a, int *lda, float *b, int *ldb, float *w, float *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
int ssygvx_(int *itype, char *jobz, char *range, char *
	uplo, int *n, float *a, int *lda, float *b, int *ldb, float *
	vl, float *vu, int *il, int *iu, float *abstol, int *m, 
	float *w, float *z__, int *ldz, float *work, int *lwork, int 
	*iwork, int *ifail, int *info);
 
int ssyrfs_(char *uplo, int *n, int *nrhs, float *a, 
	int *lda, float *af, int *ldaf, int *ipiv, float *b, 
	int *ldb, float *x, int *ldx, float *ferr, float *berr, float *
	work, int *iwork, int *info);
 
int ssysv_(char *uplo, int *n, int *nrhs, float *a, 
	int *lda, int *ipiv, float *b, int *ldb, float *work, 
	int *lwork, int *info);
 
int ssysvx_(char *fact, char *uplo, int *n, int *
	nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, 
	float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr,
	 float *berr, float *work, int *lwork, int *iwork, int *
	info);
 
int ssytd2_(char *uplo, int *n, float *a, int *lda, 
	float *d__, float *e, float *tau, int *info);
 
int ssytf2_(char *uplo, int *n, float *a, int *lda, 
	int *ipiv, int *info);
 
int ssytrd_(char *uplo, int *n, float *a, int *lda, 
	float *d__, float *e, float *tau, float *work, int *lwork, int *
	info);
 
int ssytrf_(char *uplo, int *n, float *a, int *lda, 
	int *ipiv, float *work, int *lwork, int *info);
 
int ssytri_(char *uplo, int *n, float *a, int *lda, 
	int *ipiv, float *work, int *info);
 
int ssytrs_(char *uplo, int *n, int *nrhs, float *a, 
	int *lda, int *ipiv, float *b, int *ldb, int *info);
 
int stbcon_(char *norm, char *uplo, char *diag, int *n, 
	int *kd, float *ab, int *ldab, float *rcond, float *work, 
	int *iwork, int *info);
 
int stbrfs_(char *uplo, char *trans, char *diag, int *n, 
	int *kd, int *nrhs, float *ab, int *ldab, float *b, int 
	*ldb, float *x, int *ldx, float *ferr, float *berr, float *work, 
	int *iwork, int *info);
 
int stbtrs_(char *uplo, char *trans, char *diag, int *n, 
	int *kd, int *nrhs, float *ab, int *ldab, float *b, int 
	*ldb, int *info);
 
int stgevc_(char *side, char *howmny, int *select, 
	int *n, float *a, int *lda, float *b, int *ldb, float *vl, 
	int *ldvl, float *vr, int *ldvr, int *mm, int *m, float 
	*work, int *info);
 
int stgex2_(int *wantq, int *wantz, int *n, float 
	*a, int *lda, float *b, int *ldb, float *q, int *ldq, float *
	z__, int *ldz, int *j1, int *n1, int *n2, float *work, 
	int *lwork, int *info);
 
int stgexc_(int *wantq, int *wantz, int *n, float 
	*a, int *lda, float *b, int *ldb, float *q, int *ldq, float *
	z__, int *ldz, int *ifst, int *ilst, float *work, int *
	lwork, int *info);
 
int stgsen_(int *ijob, int *wantq, int *wantz, 
	int *select, int *n, float *a, int *lda, float *b, int *
	ldb, float *alphar, float *alphai, float *beta, float *q, int *ldq, 
	float *z__, int *ldz, int *m, float *pl, float *pr, float *dif, 
	float *work, int *lwork, int *iwork, int *liwork, int *
	info);
 
int stgsja_(char *jobu, char *jobv, char *jobq, int *m, 
	int *p, int *n, int *k, int *l, float *a, int *lda,
	 float *b, int *ldb, float *tola, float *tolb, float *alpha, float *
	beta, float *u, int *ldu, float *v, int *ldv, float *q, int *
	ldq, float *work, int *ncycle, int *info);
 
int stgsna_(char *job, char *howmny, int *select, 
	int *n, float *a, int *lda, float *b, int *ldb, float *vl, 
	int *ldvl, float *vr, int *ldvr, float *s, float *dif, int *
	mm, int *m, float *work, int *lwork, int *iwork, int *
	info);
 
int stgsy2_(char *trans, int *ijob, int *m, int *
	n, float *a, int *lda, float *b, int *ldb, float *c__, int *
	ldc, float *d__, int *ldd, float *e, int *lde, float *f, int 
	*ldf, float *scale, float *rdsum, float *rdscal, int *iwork, int 
	*pq, int *info);
 
int stgsyl_(char *trans, int *ijob, int *m, int *
	n, float *a, int *lda, float *b, int *ldb, float *c__, int *
	ldc, float *d__, int *ldd, float *e, int *lde, float *f, int 
	*ldf, float *scale, float *dif, float *work, int *lwork, int *
	iwork, int *info);
 
int stpcon_(char *norm, char *uplo, char *diag, int *n, 
	float *ap, float *rcond, float *work, int *iwork, int *info);
 
int stprfs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, float *ap, float *b, int *ldb, float *x, int *ldx,
	 float *ferr, float *berr, float *work, int *iwork, int *info);
 
int stptri_(char *uplo, char *diag, int *n, float *ap, 
	int *info);
 
int stptrs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, float *ap, float *b, int *ldb, int *info);
 
int strcon_(char *norm, char *uplo, char *diag, int *n, 
	float *a, int *lda, float *rcond, float *work, int *iwork, 
	int *info);
 
int strevc_(char *side, char *howmny, int *select, 
	int *n, float *t, int *ldt, float *vl, int *ldvl, float *vr, 
	int *ldvr, int *mm, int *m, float *work, int *info);
 
int strexc_(char *compq, int *n, float *t, int *ldt, 
	float *q, int *ldq, int *ifst, int *ilst, float *work, 
	int *info);
 
int strrfs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, float *a, int *lda, float *b, int *ldb, float *x, 
	int *ldx, float *ferr, float *berr, float *work, int *iwork, 
	int *info);
 
int strsen_(char *job, char *compq, int *select, int 
	*n, float *t, int *ldt, float *q, int *ldq, float *wr, float *wi, 
	int *m, float *s, float *sep, float *work, int *lwork, int *
	iwork, int *liwork, int *info);
 
int strsna_(char *job, char *howmny, int *select, 
	int *n, float *t, int *ldt, float *vl, int *ldvl, float *vr, 
	int *ldvr, float *s, float *sep, int *mm, int *m, float *
	work, int *ldwork, int *iwork, int *info);
 
int strsyl_(char *trana, char *tranb, int *isgn, int 
	*m, int *n, float *a, int *lda, float *b, int *ldb, float *
	c__, int *ldc, float *scale, int *info);
 
int strti2_(char *uplo, char *diag, int *n, float *a, 
	int *lda, int *info);
 
int strtri_(char *uplo, char *diag, int *n, float *a, 
	int *lda, int *info);
 
int strtrs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, float *a, int *lda, float *b, int *ldb, int *
	info);
 
int stzrqf_(int *m, int *n, float *a, int *lda, 
	float *tau, int *info);
 
int stzrzf_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info);
 
int xerbla_(char *srname, int *info);
 
int zbdsqr_(char *uplo, int *n, int *ncvt, int *
	nru, int *ncc, double *d__, double *e, double_cx *vt, 
	int *ldvt, double_cx *u, int *ldu, double_cx *c__, 
	int *ldc, double *rwork, int *info);
 
int zdrot_(int *n, double_cx *cx, int *incx, 
	double_cx *cy, int *incy, double *c__, double *s);
 
int zdrscl_(int *n, double *sa, double_cx *sx, 
	int *incx);
 
int zgbbrd_(char *vect, int *m, int *n, int *ncc,
	 int *kl, int *ku, double_cx *ab, int *ldab, 
	double *d__, double *e, double_cx *q, int *ldq, 
	double_cx *pt, int *ldpt, double_cx *c__, int *ldc, 
	double_cx *work, double *rwork, int *info);
 
int zgbcon_(char *norm, int *n, int *kl, int *ku,
	 double_cx *ab, int *ldab, int *ipiv, double *anorm, 
	double *rcond, double_cx *work, double *rwork, int *
	info);
 
int zgbequ_(int *m, int *n, int *kl, int *ku,
	 double_cx *ab, int *ldab, double *r__, double *c__, 
	double *rowcnd, double *colcnd, double *amax, int *
	info);
 
int zgbrfs_(char *trans, int *n, int *kl, int *
	ku, int *nrhs, double_cx *ab, int *ldab, double_cx *
	afb, int *ldafb, int *ipiv, double_cx *b, int *ldb, 
	double_cx *x, int *ldx, double *ferr, double *berr, 
	double_cx *work, double *rwork, int *info);
 
int zgbsv_(int *n, int *kl, int *ku, int *
	nrhs, double_cx *ab, int *ldab, int *ipiv, double_cx *
	b, int *ldb, int *info);
 
int zgbsvx_(char *fact, char *trans, int *n, int *kl,
	 int *ku, int *nrhs, double_cx *ab, int *ldab, 
	double_cx *afb, int *ldafb, int *ipiv, char *equed, 
	double *r__, double *c__, double_cx *b, int *ldb, 
	double_cx *x, int *ldx, double *rcond, double *ferr, 
	double *berr, double_cx *work, double *rwork, int *
	info);
 
int zgbtf2_(int *m, int *n, int *kl, int *ku,
	 double_cx *ab, int *ldab, int *ipiv, int *info);
 
int zgbtrf_(int *m, int *n, int *kl, int *ku,
	 double_cx *ab, int *ldab, int *ipiv, int *info);
 
int zgbtrs_(char *trans, int *n, int *kl, int *
	ku, int *nrhs, double_cx *ab, int *ldab, int *ipiv, 
	double_cx *b, int *ldb, int *info);
 
int zgebak_(char *job, char *side, int *n, int *ilo, 
	int *ihi, double *scale, int *m, double_cx *v, 
	int *ldv, int *info);
 
int zgebal_(char *job, int *n, double_cx *a, int 
	*lda, int *ilo, int *ihi, double *scale, int *info);
 
int zgebd2_(int *m, int *n, double_cx *a, 
	int *lda, double *d__, double *e, double_cx *tauq, 
	double_cx *taup, double_cx *work, int *info);
 
int zgebrd_(int *m, int *n, double_cx *a, 
	int *lda, double *d__, double *e, double_cx *tauq, 
	double_cx *taup, double_cx *work, int *lwork, int *
	info);
 
int zgecon_(char *norm, int *n, double_cx *a, 
	int *lda, double *anorm, double *rcond, double_cx *
	work, double *rwork, int *info);
 
int zgeequ_(int *m, int *n, double_cx *a, 
	int *lda, double *r__, double *c__, double *rowcnd, 
	double *colcnd, double *amax, int *info);
 
int zgees_(char *jobvs, char *sort, L_fp select, int *n, 
	double_cx *a, int *lda, int *sdim, double_cx *w, 
	double_cx *vs, int *ldvs, double_cx *work, int *lwork,
	 double *rwork, int *bwork, int *info);
 
int zgeesx_(char *jobvs, char *sort, L_fp select, char *
	sense, int *n, double_cx *a, int *lda, int *sdim, 
	double_cx *w, double_cx *vs, int *ldvs, double *
	rconde, double *rcondv, double_cx *work, int *lwork, 
	double *rwork, int *bwork, int *info);
 
int zgeev_(char *jobvl, char *jobvr, int *n, 
	double_cx *a, int *lda, double_cx *w, double_cx *vl, 
	int *ldvl, double_cx *vr, int *ldvr, double_cx *work, 
	int *lwork, double *rwork, int *info);
 
int zgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, int *n, double_cx *a, int *lda, double_cx *w, 
	double_cx *vl, int *ldvl, double_cx *vr, int *ldvr, 
	int *ilo, int *ihi, double *scale, double *abnrm, 
	double *rconde, double *rcondv, double_cx *work, int *
	lwork, double *rwork, int *info);
 
int zgegs_(char *jobvsl, char *jobvsr, int *n, 
	double_cx *a, int *lda, double_cx *b, int *ldb, 
	double_cx *alpha, double_cx *beta, double_cx *vsl, 
	int *ldvsl, double_cx *vsr, int *ldvsr, double_cx *
	work, int *lwork, double *rwork, int *info);
 
int zgegv_(char *jobvl, char *jobvr, int *n, 
	double_cx *a, int *lda, double_cx *b, int *ldb, 
	double_cx *alpha, double_cx *beta, double_cx *vl, int 
	*ldvl, double_cx *vr, int *ldvr, double_cx *work, int 
	*lwork, double *rwork, int *info);
 
int zgehd2_(int *n, int *ilo, int *ihi, 
	double_cx *a, int *lda, double_cx *tau, double_cx *
	work, int *info);
 
int zgehrd_(int *n, int *ilo, int *ihi, 
	double_cx *a, int *lda, double_cx *tau, double_cx *
	work, int *lwork, int *info);
 
int zgelq2_(int *m, int *n, double_cx *a, 
	int *lda, double_cx *tau, double_cx *work, int *info);
 
int zgelqf_(int *m, int *n, double_cx *a, 
	int *lda, double_cx *tau, double_cx *work, int *lwork,
	 int *info);
 
int zgels_(char *trans, int *m, int *n, int *
	nrhs, double_cx *a, int *lda, double_cx *b, int *ldb, 
	double_cx *work, int *lwork, int *info);
 
int zgelss_(int *m, int *n, int *nrhs, 
	double_cx *a, int *lda, double_cx *b, int *ldb, double *s,
	double *rcond, int *rank, double_cx *work, int *lwork,
	double *rwork, int *info);
 
int zgelsx_(int *m, int *n, int *nrhs, 
	double_cx *a, int *lda, double_cx *b, int *ldb, 
	int *jpvt, double *rcond, int *rank, double_cx *work, 
	double *rwork, int *info);
 
int zgelsy_(int *m, int *n, int *nrhs, 
	double_cx *a, int *lda, double_cx *b, int *ldb, 
	int *jpvt, double *rcond, int *rank, double_cx *work, 
	int *lwork, double *rwork, int *info);
 
int zgeql2_(int *m, int *n, double_cx *a, 
	int *lda, double_cx *tau, double_cx *work, int *info);
 
int zgeqlf_(int *m, int *n, double_cx *a, 
	int *lda, double_cx *tau, double_cx *work, int *lwork,
	 int *info);
 
int zgeqp3_(int *m, int *n, double_cx *a, 
	int *lda, int *jpvt, double_cx *tau, double_cx *work, 
	int *lwork, double *rwork, int *info);
 
int zgeqpf_(int *m, int *n, double_cx *a, 
	int *lda, int *jpvt, double_cx *tau, double_cx *work, 
	double *rwork, int *info);
 
int zgeqr2_(int *m, int *n, double_cx *a, 
	int *lda, double_cx *tau, double_cx *work, int *info);
 
int zgeqrf_(int *m, int *n, double_cx *a, 
	int *lda, double_cx *tau, double_cx *work, int *lwork,
	 int *info);
 
int zgerfs_(char *trans, int *n, int *nrhs, 
	double_cx *a, int *lda, double_cx *af, int *ldaf, 
	int *ipiv, double_cx *b, int *ldb, double_cx *x, 
	int *ldx, double *ferr, double *berr, double_cx *work,
	 double *rwork, int *info);
 
int zgerq2_(int *m, int *n, double_cx *a, 
	int *lda, double_cx *tau, double_cx *work, int *info);
 
int zgerqf_(int *m, int *n, double_cx *a, 
	int *lda, double_cx *tau, double_cx *work, int *lwork,
	 int *info);
 
int zgesc2_(int *n, double_cx *a, int *lda, 
	double_cx *rhs, int *ipiv, int *jpiv, double *scale);

int zgesv_(int *n, int *nrhs, double_cx *a, 
	int *lda, int *ipiv, double_cx *b, int *ldb, int *
	info);
 
int zgesvd_(char *jobu, char *jobvt, int *m, int *n, 
	double_cx *a, int *lda, double *s, double_cx *u, int *
	ldu, double_cx *vt, int *ldvt, double_cx *work, int *lwork, 
	double *rwork, int *info);
 
int zgesvx_(char *fact, char *trans, int *n, int *
	nrhs, double_cx *a, int *lda, double_cx *af, int *
	ldaf, int *ipiv, char *equed, double *r__, double *c__, 
	double_cx *b, int *ldb, double_cx *x, int *ldx, 
	double *rcond, double *ferr, double *berr, double_cx *
	work, double *rwork, int *info);
 
int zgetc2_(int *n, double_cx *a, int *lda, 
	int *ipiv, int *jpiv, int *info);
 
int zgetf2_(int *m, int *n, double_cx *a, 
	int *lda, int *ipiv, int *info);
 
int zgetrf_(int *m, int *n, double_cx *a, 
	int *lda, int *ipiv, int *info);
 
int zgetri_(int *n, double_cx *a, int *lda, 
	int *ipiv, double_cx *work, int *lwork, int *info);
 
int zgetrs_(char *trans, int *n, int *nrhs, 
	double_cx *a, int *lda, int *ipiv, double_cx *b, 
	int *ldb, int *info);
 
int zggbak_(char *job, char *side, int *n, int *ilo, 
	int *ihi, double *lscale, double *rscale, int *m, 
	double_cx *v, int *ldv, int *info);
 
int zggbal_(char *job, int *n, double_cx *a, int 
	*lda, double_cx *b, int *ldb, int *ilo, int *ihi, 
	double *lscale, double *rscale, double *work, int *
	info);
 
int zgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	delctg, int *n, double_cx *a, int *lda, double_cx *b, 
	int *ldb, int *sdim, double_cx *alpha, double_cx *
	beta, double_cx *vsl, int *ldvsl, double_cx *vsr, int 
	*ldvsr, double_cx *work, int *lwork, double *rwork, 
	int *bwork, int *info);
 
int zggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	delctg, char *sense, int *n, double_cx *a, int *lda, 
	double_cx *b, int *ldb, int *sdim, double_cx *alpha, 
	double_cx *beta, double_cx *vsl, int *ldvsl, 
	double_cx *vsr, int *ldvsr, double *rconde, double *
	rcondv, double_cx *work, int *lwork, double *rwork, 
	int *iwork, int *liwork, int *bwork, int *info);
 
int zggev_(char *jobvl, char *jobvr, int *n, 
	double_cx *a, int *lda, double_cx *b, int *ldb, 
	double_cx *alpha, double_cx *beta, double_cx *vl, int 
	*ldvl, double_cx *vr, int *ldvr, double_cx *work, int 
	*lwork, double *rwork, int *info);
 
int zggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, int *n, double_cx *a, int *lda, double_cx *b, 
	int *ldb, double_cx *alpha, double_cx *beta, 
	double_cx *vl, int *ldvl, double_cx *vr, int *ldvr, 
	int *ilo, int *ihi, double *lscale, double *rscale, 
	double *abnrm, double *bbnrm, double *rconde, double *
	rcondv, double_cx *work, int *lwork, double *rwork, 
	int *iwork, int *bwork, int *info);
 
int zggglm_(int *n, int *m, int *p, 
	double_cx *a, int *lda, double_cx *b, int *ldb, 
	double_cx *d__, double_cx *x, double_cx *y, double_cx 
	*work, int *lwork, int *info);
 
int zgghrd_(char *compq, char *compz, int *n, int *
	ilo, int *ihi, double_cx *a, int *lda, double_cx *b, 
	int *ldb, double_cx *q, int *ldq, double_cx *z__, 
	int *ldz, int *info);
 
int zgglse_(int *m, int *n, int *p, 
	double_cx *a, int *lda, double_cx *b, int *ldb, 
	double_cx *c__, double_cx *d__, double_cx *x, 
	double_cx *work, int *lwork, int *info);
 
int zggqrf_(int *n, int *m, int *p, 
	double_cx *a, int *lda, double_cx *taua, double_cx *b,
	 int *ldb, double_cx *taub, double_cx *work, int *
	lwork, int *info);
 
int zggrqf_(int *m, int *p, int *n, 
	double_cx *a, int *lda, double_cx *taua, double_cx *b,
	 int *ldb, double_cx *taub, double_cx *work, int *
	lwork, int *info);
 
int zggsvd_(char *jobu, char *jobv, char *jobq, int *m, 
	int *n, int *p, int *k, int *l, double_cx *a, 
	int *lda, double_cx *b, int *ldb, double *alpha, 
	double *beta, double_cx *u, int *ldu, double_cx *v, 
	int *ldv, double_cx *q, int *ldq, double_cx *work, 
	double *rwork, int *iwork, int *info);
 
int zggsvp_(char *jobu, char *jobv, char *jobq, int *m, 
	int *p, int *n, double_cx *a, int *lda, double_cx 
	*b, int *ldb, double *tola, double *tolb, int *k, 
	int *l, double_cx *u, int *ldu, double_cx *v, int 
	*ldv, double_cx *q, int *ldq, int *iwork, double *
	rwork, double_cx *tau, double_cx *work, int *info);
 
int zgtcon_(char *norm, int *n, double_cx *dl, 
	double_cx *d__, double_cx *du, double_cx *du2, int *
	ipiv, double *anorm, double *rcond, double_cx *work, 
	int *info);
 
int zgtrfs_(char *trans, int *n, int *nrhs, 
	double_cx *dl, double_cx *d__, double_cx *du, 
	double_cx *dlf, double_cx *df, double_cx *duf, 
	double_cx *du2, int *ipiv, double_cx *b, int *ldb, 
	double_cx *x, int *ldx, double *ferr, double *berr, 
	double_cx *work, double *rwork, int *info);
 
int zgtsv_(int *n, int *nrhs, double_cx *dl, 
	double_cx *d__, double_cx *du, double_cx *b, int *ldb,
	 int *info);
 
int zgtsvx_(char *fact, char *trans, int *n, int *
	nrhs, double_cx *dl, double_cx *d__, double_cx *du, 
	double_cx *dlf, double_cx *df, double_cx *duf, 
	double_cx *du2, int *ipiv, double_cx *b, int *ldb, 
	double_cx *x, int *ldx, double *rcond, double *ferr, 
	double *berr, double_cx *work, double *rwork, int *
	info);
 
int zgttrf_(int *n, double_cx *dl, double_cx *
	d__, double_cx *du, double_cx *du2, int *ipiv, int *
	info);
 
int zgttrs_(char *trans, int *n, int *nrhs, 
	double_cx *dl, double_cx *d__, double_cx *du, 
	double_cx *du2, int *ipiv, double_cx *b, int *ldb, 
	int *info);
 
int zgtts2_(int *itrans, int *n, int *nrhs, 
	double_cx *dl, double_cx *d__, double_cx *du, 
	double_cx *du2, int *ipiv, double_cx *b, int *ldb);
 
int zhbev_(char *jobz, char *uplo, int *n, int *kd, 
	double_cx *ab, int *ldab, double *w, double_cx *z__, 
	int *ldz, double_cx *work, double *rwork, int *info);
 
int zhbevd_(char *jobz, char *uplo, int *n, int *kd, 
	double_cx *ab, int *ldab, double *w, double_cx *z__, 
	int *ldz, double_cx *work, int *lwork, double *rwork, 
	int *lrwork, int *iwork, int *liwork, int *info);
 
int zhbevx_(char *jobz, char *range, char *uplo, int *n, 
	int *kd, double_cx *ab, int *ldab, double_cx *q, 
	int *ldq, double *vl, double *vu, int *il, int *
	iu, double *abstol, int *m, double *w, double_cx *z__,
	 int *ldz, double_cx *work, double *rwork, int *iwork,
	 int *ifail, int *info);
 
int zhbgst_(char *vect, char *uplo, int *n, int *ka, 
	int *kb, double_cx *ab, int *ldab, double_cx *bb, 
	int *ldbb, double_cx *x, int *ldx, double_cx *work, 
	double *rwork, int *info);
 
int zhbgv_(char *jobz, char *uplo, int *n, int *ka, 
	int *kb, double_cx *ab, int *ldab, double_cx *bb, 
	int *ldbb, double *w, double_cx *z__, int *ldz, 
	double_cx *work, double *rwork, int *info);
 
int zhbgvx_(char *jobz, char *range, char *uplo, int *n, 
	int *ka, int *kb, double_cx *ab, int *ldab, 
	double_cx *bb, int *ldbb, double_cx *q, int *ldq, 
	double *vl, double *vu, int *il, int *iu, double *
	abstol, int *m, double *w, double_cx *z__, int *ldz, 
	double_cx *work, double *rwork, int *iwork, int *
	ifail, int *info);
 
int zhbtrd_(char *vect, char *uplo, int *n, int *kd, 
	double_cx *ab, int *ldab, double *d__, double *e, 
	double_cx *q, int *ldq, double_cx *work, int *info);
 
int zhecon_(char *uplo, int *n, double_cx *a, 
	int *lda, int *ipiv, double *anorm, double *rcond, 
	double_cx *work, int *info);
 
int zheev_(char *jobz, char *uplo, int *n, double_cx 
	*a, int *lda, double *w, double_cx *work, int *lwork, 
	double *rwork, int *info);
 
int zheevd_(char *jobz, char *uplo, int *n, 
	double_cx *a, int *lda, double *w, double_cx *work, 
	int *lwork, double *rwork, int *lrwork, int *iwork, 
	int *liwork, int *info);
 
int zheevr_(char *jobz, char *range, char *uplo, int *n, 
	double_cx *a, int *lda, double *vl, double *vu, 
	int *il, int *iu, double *abstol, int *m, double *
	w, double_cx *z__, int *ldz, int *isuppz, double_cx *
	work, int *lwork, double *rwork, int *lrwork, int *
	iwork, int *liwork, int *info);
 
int zheevx_(char *jobz, char *range, char *uplo, int *n, 
	double_cx *a, int *lda, double *vl, double *vu, 
	int *il, int *iu, double *abstol, int *m, double *
	w, double_cx *z__, int *ldz, double_cx *work, int *
	lwork, double *rwork, int *iwork, int *ifail, int *
	info);
 
int zhegs2_(int *itype, char *uplo, int *n, 
	double_cx *a, int *lda, double_cx *b, int *ldb, 
	int *info);
 
int zhegst_(int *itype, char *uplo, int *n, 
	double_cx *a, int *lda, double_cx *b, int *ldb, 
	int *info);
 
int zhegv_(int *itype, char *jobz, char *uplo, int *
	n, double_cx *a, int *lda, double_cx *b, int *ldb, 
	double *w, double_cx *work, int *lwork, double *rwork,
	 int *info);
 
int zhegvd_(int *itype, char *jobz, char *uplo, int *
	n, double_cx *a, int *lda, double_cx *b, int *ldb, 
	double *w, double_cx *work, int *lwork, double *rwork,
	 int *lrwork, int *iwork, int *liwork, int *info);
 
int zhegvx_(int *itype, char *jobz, char *range, char *
	uplo, int *n, double_cx *a, int *lda, double_cx *b, 
	int *ldb, double *vl, double *vu, int *il, int *
	iu, double *abstol, int *m, double *w, double_cx *z__,
	 int *ldz, double_cx *work, int *lwork, double *rwork,
	 int *iwork, int *ifail, int *info);
 
int zherfs_(char *uplo, int *n, int *nrhs, 
	double_cx *a, int *lda, double_cx *af, int *ldaf, 
	int *ipiv, double_cx *b, int *ldb, double_cx *x, 
	int *ldx, double *ferr, double *berr, double_cx *work,
	 double *rwork, int *info);
 
int zhesv_(char *uplo, int *n, int *nrhs, 
	double_cx *a, int *lda, int *ipiv, double_cx *b, 
	int *ldb, double_cx *work, int *lwork, int *info);
 
int zhesvx_(char *fact, char *uplo, int *n, int *
	nrhs, double_cx *a, int *lda, double_cx *af, int *
	ldaf, int *ipiv, double_cx *b, int *ldb, double_cx *x,
	 int *ldx, double *rcond, double *ferr, double *berr, 
	double_cx *work, int *lwork, double *rwork, int *info);
 
int zhetf2_(char *uplo, int *n, double_cx *a, 
	int *lda, int *ipiv, int *info);
 
int zhetrd_(char *uplo, int *n, double_cx *a, 
	int *lda, double *d__, double *e, double_cx *tau, 
	double_cx *work, int *lwork, int *info);
 
int zhetrf_(char *uplo, int *n, double_cx *a, 
	int *lda, int *ipiv, double_cx *work, int *lwork, 
	int *info);
 
int zhetri_(char *uplo, int *n, double_cx *a, 
	int *lda, int *ipiv, double_cx *work, int *info);
 
int zhetrs_(char *uplo, int *n, int *nrhs, 
	double_cx *a, int *lda, int *ipiv, double_cx *b, 
	int *ldb, int *info);
 
int zhgeqz_(char *job, char *compq, char *compz, int *n, 
	int *ilo, int *ihi, double_cx *a, int *lda, 
	double_cx *b, int *ldb, double_cx *alpha, double_cx *
	beta, double_cx *q, int *ldq, double_cx *z__, int *
	ldz, double_cx *work, int *lwork, double *rwork, int *
	info);
 
int zhpcon_(char *uplo, int *n, double_cx *ap, 
	int *ipiv, double *anorm, double *rcond, double_cx *
	work, int *info);
 
int zhpev_(char *jobz, char *uplo, int *n, double_cx 
	*ap, double *w, double_cx *z__, int *ldz, double_cx *
	work, double *rwork, int *info);
 
int zhpevd_(char *jobz, char *uplo, int *n, 
	double_cx *ap, double *w, double_cx *z__, int *ldz, 
	double_cx *work, int *lwork, double *rwork, int *
	lrwork, int *iwork, int *liwork, int *info);
 
int zhpevx_(char *jobz, char *range, char *uplo, int *n, 
	double_cx *ap, double *vl, double *vu, int *il, 
	int *iu, double *abstol, int *m, double *w, 
	double_cx *z__, int *ldz, double_cx *work, double *
	rwork, int *iwork, int *ifail, int *info);
 
int zhpgst_(int *itype, char *uplo, int *n, 
	double_cx *ap, double_cx *bp, int *info);
 
int zhpgv_(int *itype, char *jobz, char *uplo, int *
	n, double_cx *ap, double_cx *bp, double *w, double_cx 
	*z__, int *ldz, double_cx *work, double *rwork, int *
	info);
 
int zhpgvd_(int *itype, char *jobz, char *uplo, int *
	n, double_cx *ap, double_cx *bp, double *w, double_cx 
	*z__, int *ldz, double_cx *work, int *lwork, double *
	rwork, int *lrwork, int *iwork, int *liwork, int *
	info);
 
int zhpgvx_(int *itype, char *jobz, char *range, char *
	uplo, int *n, double_cx *ap, double_cx *bp, double *
	vl, double *vu, int *il, int *iu, double *abstol, 
	int *m, double *w, double_cx *z__, int *ldz, 
	double_cx *work, double *rwork, int *iwork, int *
	ifail, int *info);
 
int zhprfs_(char *uplo, int *n, int *nrhs, 
	double_cx *ap, double_cx *afp, int *ipiv, double_cx *
	b, int *ldb, double_cx *x, int *ldx, double *ferr, 
	double *berr, double_cx *work, double *rwork, int *
	info);
 
int zhpsv_(char *uplo, int *n, int *nrhs, 
	double_cx *ap, int *ipiv, double_cx *b, int *ldb, 
	int *info);
 
int zhpsvx_(char *fact, char *uplo, int *n, int *
	nrhs, double_cx *ap, double_cx *afp, int *ipiv, 
	double_cx *b, int *ldb, double_cx *x, int *ldx, 
	double *rcond, double *ferr, double *berr, double_cx *
	work, double *rwork, int *info);
 
int zhptrd_(char *uplo, int *n, double_cx *ap, 
	double *d__, double *e, double_cx *tau, int *info);
 
int zhptrf_(char *uplo, int *n, double_cx *ap, 
	int *ipiv, int *info);
 
int zhptri_(char *uplo, int *n, double_cx *ap, 
	int *ipiv, double_cx *work, int *info);
 
int zhptrs_(char *uplo, int *n, int *nrhs, 
	double_cx *ap, int *ipiv, double_cx *b, int *ldb, 
	int *info);
 
int zhsein_(char *side, char *eigsrc, char *initv, int *
	select, int *n, double_cx *h__, int *ldh, double_cx *
	w, double_cx *vl, int *ldvl, double_cx *vr, int *ldvr,
	 int *mm, int *m, double_cx *work, double *rwork, 
	int *ifaill, int *ifailr, int *info);
 
int zhseqr_(char *job, char *compz, int *n, int *ilo,
	 int *ihi, double_cx *h__, int *ldh, double_cx *w, 
	double_cx *z__, int *ldz, double_cx *work, int *lwork,
	 int *info);
 
int zlabrd_(int *m, int *n, int *nb, 
	double_cx *a, int *lda, double *d__, double *e, 
	double_cx *tauq, double_cx *taup, double_cx *x, int *
	ldx, double_cx *y, int *ldy);
 
int zlacgv_(int *n, double_cx *x, int *incx);
 
int zlacon_(int *n, double_cx *v, double_cx *x, 
	double *est, int *kase);
 
int zlacp2_(char *uplo, int *m, int *n, double *
	a, int *lda, double_cx *b, int *ldb);
 
int zlacpy_(char *uplo, int *m, int *n, 
	double_cx *a, int *lda, double_cx *b, int *ldb);
 
int zlacrm_(int *m, int *n, double_cx *a, 
	int *lda, double *b, int *ldb, double_cx *c__, 
	int *ldc, double *rwork);
 
int zlacrt_(int *n, double_cx *cx, int *incx, 
	double_cx *cy, int *incy, double_cx *c__, double_cx *
	s);
 
int zlaed0_(int *qsiz, int *n, double *d__, 
	double *e, double_cx *q, int *ldq, double_cx *qstore, 
	int *ldqs, double *rwork, int *iwork, int *info);
 
int zlaed7_(int *n, int *cutpnt, int *qsiz, 
	int *tlvls, int *curlvl, int *curpbm, double *d__, 
	double_cx *q, int *ldq, double *rho, int *indxq, 
	double *qstore, int *qptr, int *prmptr, int *perm, 
	int *givptr, int *givcol, double *givnum, double_cx *
	work, double *rwork, int *iwork, int *info);
 
int zlaed8_(int *k, int *n, int *qsiz, 
	double_cx *q, int *ldq, double *d__, double *rho, 
	int *cutpnt, double *z__, double *dlamda, double_cx *
	q2, int *ldq2, double *w, int *indxp, int *indx, 
	int *indxq, int *perm, int *givptr, int *givcol, 
	double *givnum, int *info);
 
int zlaein_(int *rightv, int *noinit, int *n, 
	double_cx *h__, int *ldh, double_cx *w, double_cx *v, 
	double_cx *b, int *ldb, double *rwork, double *eps3, 
	double *smlnum, int *info);
 
int zlaesy_(double_cx *a, double_cx *b, 
	double_cx *c__, double_cx *rt1, double_cx *rt2, 
	double_cx *evscal, double_cx *cs1, double_cx *sn1);
 
int zlaev2_(double_cx *a, double_cx *b, 
	double_cx *c__, double *rt1, double *rt2, double *cs1,
	 double_cx *sn1);
 
int zlags2_(int *upper, double *a1, double_cx *
	a2, double *a3, double *b1, double_cx *b2, double *b3,
	 double *csu, double_cx *snu, double *csv, double_cx *
	snv, double *csq, double_cx *snq);
 
int zlagtm_(char *trans, int *n, int *nrhs, 
	double *alpha, double_cx *dl, double_cx *d__, 
	double_cx *du, double_cx *x, int *ldx, double *beta, 
	double_cx *b, int *ldb);
 
int zlahef_(char *uplo, int *n, int *nb, int *kb,
	 double_cx *a, int *lda, int *ipiv, double_cx *w, 
	int *ldw, int *info);
 
int zlahqr_(int *wantt, int *wantz, int *n, 
	int *ilo, int *ihi, double_cx *h__, int *ldh, 
	double_cx *w, int *iloz, int *ihiz, double_cx *z__, 
	int *ldz, int *info);
 
int zlahrd_(int *n, int *k, int *nb, 
	double_cx *a, int *lda, double_cx *tau, double_cx *t, 
	int *ldt, double_cx *y, int *ldy);
 
int zlaic1_(int *job, int *j, double_cx *x, 
	double *sest, double_cx *w, double_cx *gamma, double *
	sestpr, double_cx *s, double_cx *c__);
 
int zlals0_(int *icompq, int *nl, int *nr, 
	int *sqre, int *nrhs, double_cx *b, int *ldb, 
	double_cx *bx, int *ldbx, int *perm, int *givptr, 
	int *givcol, int *ldgcol, double *givnum, int *ldgnum,
	 double *poles, double *difl, double *difr, double *
	z__, int *k, double *c__, double *s, double *rwork, 
	int *info);
 
int zlalsa_(int *icompq, int *smlsiz, int *n, 
	int *nrhs, double_cx *b, int *ldb, double_cx *bx, 
	int *ldbx, double *u, int *ldu, double *vt, int *
	k, double *difl, double *difr, double *z__, double *
	poles, int *givptr, int *givcol, int *ldgcol, int *
	perm, double *givnum, double *c__, double *s, double *
	rwork, int *iwork, int *info);
 
int zlapll_(int *n, double_cx *x, int *incx, 
	double_cx *y, int *incy, double *ssmin);
 
int zlapmt_(int *forwrd, int *m, int *n, 
	double_cx *x, int *ldx, int *k);
 
int zlaqgb_(int *m, int *n, int *kl, int *ku,
	 double_cx *ab, int *ldab, double *r__, double *c__, 
	double *rowcnd, double *colcnd, double *amax, char *equed);
 
int zlaqge_(int *m, int *n, double_cx *a, 
	int *lda, double *r__, double *c__, double *rowcnd, 
	double *colcnd, double *amax, char *equed);
 
int zlaqhb_(char *uplo, int *n, int *kd, 
	double_cx *ab, int *ldab, double *s, double *scond, 
	double *amax, char *equed);
 
int zlaqhe_(char *uplo, int *n, double_cx *a, 
	int *lda, double *s, double *scond, double *amax, 
	char *equed);
 
int zlaqhp_(char *uplo, int *n, double_cx *ap, 
	double *s, double *scond, double *amax, char *equed);
 
int zlaqp2_(int *m, int *n, int *offset, 
	double_cx *a, int *lda, int *jpvt, double_cx *tau, 
	double *vn1, double *vn2, double_cx *work);
 
int zlaqps_(int *m, int *n, int *offset, int 
	*nb, int *kb, double_cx *a, int *lda, int *jpvt, 
	double_cx *tau, double *vn1, double *vn2, double_cx *
	auxv, double_cx *f, int *ldf);
 
int zlaqsb_(char *uplo, int *n, int *kd, 
	double_cx *ab, int *ldab, double *s, double *scond, 
	double *amax, char *equed);
 
int zlaqsp_(char *uplo, int *n, double_cx *ap, 
	double *s, double *scond, double *amax, char *equed);
 
int zlaqsy_(char *uplo, int *n, double_cx *a, 
	int *lda, double *s, double *scond, double *amax, 
	char *equed);
 
int zlar1v_(int *n, int *b1, int *bn, double 
	*sigma, double *d__, double *l, double *ld, double *
	lld, double *gersch, double_cx *z__, double *ztz, 
	double *mingma, int *r__, int *isuppz, double *work);
 
int zlar2v_(int *n, double_cx *x, double_cx *y, 
	double_cx *z__, int *incx, double *c__, double_cx *s, 
	int *incc);
 
int zlarcm_(int *m, int *n, double *a, int *
	lda, double_cx *b, int *ldb, double_cx *c__, int *ldc,
	 double *rwork);
 
int zlarf_(char *side, int *m, int *n, double_cx 
	*v, int *incv, double_cx *tau, double_cx *c__, int *
	ldc, double_cx *work);
 
int zlarfb_(char *side, char *trans, char *direct, char *
	storev, int *m, int *n, int *k, double_cx *v, int 
	*ldv, double_cx *t, int *ldt, double_cx *c__, int *
	ldc, double_cx *work, int *ldwork);
 
int zlarfg_(int *n, double_cx *alpha, double_cx *
	x, int *incx, double_cx *tau);
 
int zlarft_(char *direct, char *storev, int *n, int *
	k, double_cx *v, int *ldv, double_cx *tau, double_cx *
	t, int *ldt);
 
int zlarfx_(char *side, int *m, int *n, 
	double_cx *v, double_cx *tau, double_cx *c__, int *
	ldc, double_cx *work);
 
int zlargv_(int *n, double_cx *x, int *incx, 
	double_cx *y, int *incy, double *c__, int *incc);
 
int zlarnv_(int *idist, int *iseed, int *n, 
	double_cx *x);
 
int zlarrv_(int *n, double *d__, double *l, 
	int *isplit, int *m, double *w, int *iblock, 
	double *gersch, double *tol, double_cx *z__, int *ldz,
	 int *isuppz, double *work, int *iwork, int *info);
 
int zlartg_(double_cx *f, double_cx *g, double *
	cs, double_cx *sn, double_cx *r__);
 
int zlartv_(int *n, double_cx *x, int *incx, 
	double_cx *y, int *incy, double *c__, double_cx *s, 
	int *incc);
 
int zlarz_(char *side, int *m, int *n, int *l, 
	double_cx *v, int *incv, double_cx *tau, double_cx *
	c__, int *ldc, double_cx *work);
 
int zlarzb_(char *side, char *trans, char *direct, char *
	storev, int *m, int *n, int *k, int *l, double_cx 
	*v, int *ldv, double_cx *t, int *ldt, double_cx *c__, 
	int *ldc, double_cx *work, int *ldwork);
 
int zlarzt_(char *direct, char *storev, int *n, int *
	k, double_cx *v, int *ldv, double_cx *tau, double_cx *
	t, int *ldt);
 
int zlascl_(char *type__, int *kl, int *ku, 
	double *cfrom, double *cto, int *m, int *n, 
	double_cx *a, int *lda, int *info);
 
int zlaset_(char *uplo, int *m, int *n, 
	double_cx *alpha, double_cx *beta, double_cx *a, int *
	lda);
 
int zlasr_(char *side, char *pivot, char *direct, int *m,
	 int *n, double *c__, double *s, double_cx *a, 
	int *lda);
 
int zlassq_(int *n, double_cx *x, int *incx, 
	double *scale, double *sumsq);
 
int zlaswp_(int *n, double_cx *a, int *lda, 
	int *k1, int *k2, int *ipiv, int *incx);
 
int zlasyf_(char *uplo, int *n, int *nb, int *kb,
	 double_cx *a, int *lda, int *ipiv, double_cx *w, 
	int *ldw, int *info);
 
int zlatbs_(char *uplo, char *trans, char *diag, char *
	normin, int *n, int *kd, double_cx *ab, int *ldab, 
	double_cx *x, double *scale, double *cnorm, int *info);
 
int zlatdf_(int *ijob, int *n, double_cx *z__, 
	int *ldz, double_cx *rhs, double *rdsum, double *
	rdscal, int *ipiv, int *jpiv);
 
int zlatps_(char *uplo, char *trans, char *diag, char *
	normin, int *n, double_cx *ap, double_cx *x, double *
	scale, double *cnorm, int *info);
 
int zlatrd_(char *uplo, int *n, int *nb, 
	double_cx *a, int *lda, double *e, double_cx *tau, 
	double_cx *w, int *ldw);
 
int zlatrs_(char *uplo, char *trans, char *diag, char *
	normin, int *n, double_cx *a, int *lda, double_cx *x, 
	double *scale, double *cnorm, int *info);
 
int zlatrz_(int *m, int *n, int *l, 
	double_cx *a, int *lda, double_cx *tau, double_cx *
	work);
 
int zlatzm_(char *side, int *m, int *n, 
	double_cx *v, int *incv, double_cx *tau, double_cx *
	c1, double_cx *c2, int *ldc, double_cx *work);
 
int zlauu2_(char *uplo, int *n, double_cx *a, 
	int *lda, int *info);
 
int zlauum_(char *uplo, int *n, double_cx *a, 
	int *lda, int *info);
 
int zpbcon_(char *uplo, int *n, int *kd, 
	double_cx *ab, int *ldab, double *anorm, double *
	rcond, double_cx *work, double *rwork, int *info);
 
int zpbequ_(char *uplo, int *n, int *kd, 
	double_cx *ab, int *ldab, double *s, double *scond, 
	double *amax, int *info);
 
int zpbrfs_(char *uplo, int *n, int *kd, int *
	nrhs, double_cx *ab, int *ldab, double_cx *afb, int *
	ldafb, double_cx *b, int *ldb, double_cx *x, int *ldx,
	 double *ferr, double *berr, double_cx *work, double *
	rwork, int *info);
 
int zpbstf_(char *uplo, int *n, int *kd, 
	double_cx *ab, int *ldab, int *info);
 
int zpbsv_(char *uplo, int *n, int *kd, int *
	nrhs, double_cx *ab, int *ldab, double_cx *b, int *
	ldb, int *info);
 
int zpbsvx_(char *fact, char *uplo, int *n, int *kd, 
	int *nrhs, double_cx *ab, int *ldab, double_cx *afb, 
	int *ldafb, char *equed, double *s, double_cx *b, int 
	*ldb, double_cx *x, int *ldx, double *rcond, double *
	ferr, double *berr, double_cx *work, double *rwork, 
	int *info);
 
int zpbtf2_(char *uplo, int *n, int *kd, 
	double_cx *ab, int *ldab, int *info);
 
int zpbtrf_(char *uplo, int *n, int *kd, 
	double_cx *ab, int *ldab, int *info);
 
int zpbtrs_(char *uplo, int *n, int *kd, int *
	nrhs, double_cx *ab, int *ldab, double_cx *b, int *
	ldb, int *info);
 
int zpocon_(char *uplo, int *n, double_cx *a, 
	int *lda, double *anorm, double *rcond, double_cx *
	work, double *rwork, int *info);
 
int zpoequ_(int *n, double_cx *a, int *lda, 
	double *s, double *scond, double *amax, int *info);
 
int zporfs_(char *uplo, int *n, int *nrhs, 
	double_cx *a, int *lda, double_cx *af, int *ldaf, 
	double_cx *b, int *ldb, double_cx *x, int *ldx, 
	double *ferr, double *berr, double_cx *work, double *
	rwork, int *info);
 
int zposv_(char *uplo, int *n, int *nrhs, 
	double_cx *a, int *lda, double_cx *b, int *ldb, 
	int *info);
 
int zposvx_(char *fact, char *uplo, int *n, int *
	nrhs, double_cx *a, int *lda, double_cx *af, int *
	ldaf, char *equed, double *s, double_cx *b, int *ldb, 
	double_cx *x, int *ldx, double *rcond, double *ferr, 
	double *berr, double_cx *work, double *rwork, int *
	info);
 
int zpotf2_(char *uplo, int *n, double_cx *a, 
	int *lda, int *info);
 
int zpotrf_(char *uplo, int *n, double_cx *a, 
	int *lda, int *info);
 
int zpotri_(char *uplo, int *n, double_cx *a, 
	int *lda, int *info);
 
int zpotrs_(char *uplo, int *n, int *nrhs, 
	double_cx *a, int *lda, double_cx *b, int *ldb, 
	int *info);
 
int zppcon_(char *uplo, int *n, double_cx *ap, 
	double *anorm, double *rcond, double_cx *work, double 
	*rwork, int *info);
 
int zppequ_(char *uplo, int *n, double_cx *ap, 
	double *s, double *scond, double *amax, int *info);
 
int zpprfs_(char *uplo, int *n, int *nrhs, 
	double_cx *ap, double_cx *afp, double_cx *b, int *ldb,
	 double_cx *x, int *ldx, double *ferr, double *berr, 
	double_cx *work, double *rwork, int *info);
 
int zppsv_(char *uplo, int *n, int *nrhs, 
	double_cx *ap, double_cx *b, int *ldb, int *info);
 
int zppsvx_(char *fact, char *uplo, int *n, int *
	nrhs, double_cx *ap, double_cx *afp, char *equed, double *
	s, double_cx *b, int *ldb, double_cx *x, int *ldx, 
	double *rcond, double *ferr, double *berr, double_cx *
	work, double *rwork, int *info);
 
int zpptrf_(char *uplo, int *n, double_cx *ap, 
	int *info);
 
int zpptri_(char *uplo, int *n, double_cx *ap, 
	int *info);
 
int zpptrs_(char *uplo, int *n, int *nrhs, 
	double_cx *ap, double_cx *b, int *ldb, int *info);
 
int zptcon_(int *n, double *d__, double_cx *e, 
	double *anorm, double *rcond, double *rwork, int *
	info);
 
int zptrfs_(char *uplo, int *n, int *nrhs, 
	double *d__, double_cx *e, double *df, double_cx *ef, 
	double_cx *b, int *ldb, double_cx *x, int *ldx, 
	double *ferr, double *berr, double_cx *work, double *
	rwork, int *info);
 
int zptsv_(int *n, int *nrhs, double *d__, 
	double_cx *e, double_cx *b, int *ldb, int *info);
 
int zptsvx_(char *fact, int *n, int *nrhs, 
	double *d__, double_cx *e, double *df, double_cx *ef, 
	double_cx *b, int *ldb, double_cx *x, int *ldx, 
	double *rcond, double *ferr, double *berr, double_cx *
	work, double *rwork, int *info);
 
int zpttrf_(int *n, double *d__, double_cx *e, 
	int *info);
 
int zpttrs_(char *uplo, int *n, int *nrhs, 
	double *d__, double_cx *e, double_cx *b, int *ldb, 
	int *info);
 
int zptts2_(int *iuplo, int *n, int *nrhs, 
	double *d__, double_cx *e, double_cx *b, int *ldb);
 
int zrot_(int *n, double_cx *cx, int *incx, 
	double_cx *cy, int *incy, double *c__, double_cx *s);
 
int zspcon_(char *uplo, int *n, double_cx *ap, 
	int *ipiv, double *anorm, double *rcond, double_cx *
	work, int *info);
 
int zspmv_(char *uplo, int *n, double_cx *alpha, 
	double_cx *ap, double_cx *x, int *incx, double_cx *
	beta, double_cx *y, int *incy);
 
int zspr_(char *uplo, int *n, double_cx *alpha, 
	double_cx *x, int *incx, double_cx *ap);
 
int zsprfs_(char *uplo, int *n, int *nrhs, 
	double_cx *ap, double_cx *afp, int *ipiv, double_cx *
	b, int *ldb, double_cx *x, int *ldx, double *ferr, 
	double *berr, double_cx *work, double *rwork, int *
	info);
 
int zspsv_(char *uplo, int *n, int *nrhs, 
	double_cx *ap, int *ipiv, double_cx *b, int *ldb, 
	int *info);
 
int zspsvx_(char *fact, char *uplo, int *n, int *
	nrhs, double_cx *ap, double_cx *afp, int *ipiv, 
	double_cx *b, int *ldb, double_cx *x, int *ldx, 
	double *rcond, double *ferr, double *berr, double_cx *
	work, double *rwork, int *info);
 
int zsptrf_(char *uplo, int *n, double_cx *ap, 
	int *ipiv, int *info);
 
int zsptri_(char *uplo, int *n, double_cx *ap, 
	int *ipiv, double_cx *work, int *info);
 
int zsptrs_(char *uplo, int *n, int *nrhs, 
	double_cx *ap, int *ipiv, double_cx *b, int *ldb, 
	int *info);
 
int zstedc_(char *compz, int *n, double *d__, 
	double *e, double_cx *z__, int *ldz, double_cx *work, 
	int *lwork, double *rwork, int *lrwork, int *iwork, 
	int *liwork, int *info);
 
int zstein_(int *n, double *d__, double *e, 
	int *m, double *w, int *iblock, int *isplit, 
	double_cx *z__, int *ldz, double *work, int *iwork, 
	int *ifail, int *info);
 
int zsteqr_(char *compz, int *n, double *d__, 
	double *e, double_cx *z__, int *ldz, double *work, 
	int *info);
 
int zsycon_(char *uplo, int *n, double_cx *a, 
	int *lda, int *ipiv, double *anorm, double *rcond, 
	double_cx *work, int *info);
 
int zsymv_(char *uplo, int *n, double_cx *alpha, 
	double_cx *a, int *lda, double_cx *x, int *incx, 
	double_cx *beta, double_cx *y, int *incy);
 
int zsyr_(char *uplo, int *n, double_cx *alpha, 
	double_cx *x, int *incx, double_cx *a, int *lda);
 
int zsyrfs_(char *uplo, int *n, int *nrhs, 
	double_cx *a, int *lda, double_cx *af, int *ldaf, 
	int *ipiv, double_cx *b, int *ldb, double_cx *x, 
	int *ldx, double *ferr, double *berr, double_cx *work,
	 double *rwork, int *info);
 
int zsysv_(char *uplo, int *n, int *nrhs, 
	double_cx *a, int *lda, int *ipiv, double_cx *b, 
	int *ldb, double_cx *work, int *lwork, int *info);
 
int zsysvx_(char *fact, char *uplo, int *n, int *
	nrhs, double_cx *a, int *lda, double_cx *af, int *
	ldaf, int *ipiv, double_cx *b, int *ldb, double_cx *x,
	 int *ldx, double *rcond, double *ferr, double *berr, 
	double_cx *work, int *lwork, double *rwork, int *info);
 
int zsytf2_(char *uplo, int *n, double_cx *a, 
	int *lda, int *ipiv, int *info);
 
int zsytrf_(char *uplo, int *n, double_cx *a, 
	int *lda, int *ipiv, double_cx *work, int *lwork, 
	int *info);
 
int zsytri_(char *uplo, int *n, double_cx *a, 
	int *lda, int *ipiv, double_cx *work, int *info);
 
int zsytrs_(char *uplo, int *n, int *nrhs, 
	double_cx *a, int *lda, int *ipiv, double_cx *b, 
	int *ldb, int *info);
 
int ztbcon_(char *norm, char *uplo, char *diag, int *n, 
	int *kd, double_cx *ab, int *ldab, double *rcond, 
	double_cx *work, double *rwork, int *info);
 
int ztbrfs_(char *uplo, char *trans, char *diag, int *n, 
	int *kd, int *nrhs, double_cx *ab, int *ldab, 
	double_cx *b, int *ldb, double_cx *x, int *ldx, 
	double *ferr, double *berr, double_cx *work, double *
	rwork, int *info);
 
int ztbtrs_(char *uplo, char *trans, char *diag, int *n, 
	int *kd, int *nrhs, double_cx *ab, int *ldab, 
	double_cx *b, int *ldb, int *info);
 
int ztgevc_(char *side, char *howmny, int *select, 
	int *n, double_cx *a, int *lda, double_cx *b, int 
	*ldb, double_cx *vl, int *ldvl, double_cx *vr, int *
	ldvr, int *mm, int *m, double_cx *work, double *rwork,
	 int *info);
 
int ztgex2_(int *wantq, int *wantz, int *n, 
	double_cx *a, int *lda, double_cx *b, int *ldb, 
	double_cx *q, int *ldq, double_cx *z__, int *ldz, 
	int *j1, int *info);
 
int ztgexc_(int *wantq, int *wantz, int *n, 
	double_cx *a, int *lda, double_cx *b, int *ldb, 
	double_cx *q, int *ldq, double_cx *z__, int *ldz, 
	int *ifst, int *ilst, int *info);
 
int ztgsen_(int *ijob, int *wantq, int *wantz, 
	int *select, int *n, double_cx *a, int *lda, 
	double_cx *b, int *ldb, double_cx *alpha, double_cx *
	beta, double_cx *q, int *ldq, double_cx *z__, int *
	ldz, int *m, double *pl, double *pr, double *dif, 
	double_cx *work, int *lwork, int *iwork, int *liwork, 
	int *info);
 
int ztgsja_(char *jobu, char *jobv, char *jobq, int *m, 
	int *p, int *n, int *k, int *l, double_cx *a, 
	int *lda, double_cx *b, int *ldb, double *tola, 
	double *tolb, double *alpha, double *beta, double_cx *
	u, int *ldu, double_cx *v, int *ldv, double_cx *q, 
	int *ldq, double_cx *work, int *ncycle, int *info);
 
int ztgsna_(char *job, char *howmny, int *select, 
	int *n, double_cx *a, int *lda, double_cx *b, int 
	*ldb, double_cx *vl, int *ldvl, double_cx *vr, int *
	ldvr, double *s, double *dif, int *mm, int *m, 
	double_cx *work, int *lwork, int *iwork, int *info);
 
int ztgsy2_(char *trans, int *ijob, int *m, int *
	n, double_cx *a, int *lda, double_cx *b, int *ldb, 
	double_cx *c__, int *ldc, double_cx *d__, int *ldd, 
	double_cx *e, int *lde, double_cx *f, int *ldf, 
	double *scale, double *rdsum, double *rdscal, int *
	info);
 
int ztgsyl_(char *trans, int *ijob, int *m, int *
	n, double_cx *a, int *lda, double_cx *b, int *ldb, 
	double_cx *c__, int *ldc, double_cx *d__, int *ldd, 
	double_cx *e, int *lde, double_cx *f, int *ldf, 
	double *scale, double *dif, double_cx *work, int *
	lwork, int *iwork, int *info);
 
int ztpcon_(char *norm, char *uplo, char *diag, int *n, 
	double_cx *ap, double *rcond, double_cx *work, double 
	*rwork, int *info);
 
int ztprfs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, double_cx *ap, double_cx *b, int *ldb, 
	double_cx *x, int *ldx, double *ferr, double *berr, 
	double_cx *work, double *rwork, int *info);
 
int ztptri_(char *uplo, char *diag, int *n, 
	double_cx *ap, int *info);
 
int ztptrs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, double_cx *ap, double_cx *b, int *ldb, 
	int *info);
 
int ztrcon_(char *norm, char *uplo, char *diag, int *n, 
	double_cx *a, int *lda, double *rcond, double_cx *
	work, double *rwork, int *info);
 
int ztrevc_(char *side, char *howmny, int *select, 
	int *n, double_cx *t, int *ldt, double_cx *vl, 
	int *ldvl, double_cx *vr, int *ldvr, int *mm, int 
	*m, double_cx *work, double *rwork, int *info);
 
int ztrexc_(char *compq, int *n, double_cx *t, 
	int *ldt, double_cx *q, int *ldq, int *ifst, int *
	ilst, int *info);
 
int ztrrfs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, double_cx *a, int *lda, double_cx *b, 
	int *ldb, double_cx *x, int *ldx, double *ferr, 
	double *berr, double_cx *work, double *rwork, int *
	info);
 
int ztrsen_(char *job, char *compq, int *select, int 
	*n, double_cx *t, int *ldt, double_cx *q, int *ldq, 
	double_cx *w, int *m, double *s, double *sep, 
	double_cx *work, int *lwork, int *info);
 
int ztrsna_(char *job, char *howmny, int *select, 
	int *n, double_cx *t, int *ldt, double_cx *vl, 
	int *ldvl, double_cx *vr, int *ldvr, double *s, 
	double *sep, int *mm, int *m, double_cx *work, 
	int *ldwork, double *rwork, int *info);
 
int ztrsyl_(char *trana, char *tranb, int *isgn, int 
	*m, int *n, double_cx *a, int *lda, double_cx *b, 
	int *ldb, double_cx *c__, int *ldc, double *scale, 
	int *info);
 
int ztrti2_(char *uplo, char *diag, int *n, 
	double_cx *a, int *lda, int *info);
 
int ztrtri_(char *uplo, char *diag, int *n, 
	double_cx *a, int *lda, int *info);
 
int ztrtrs_(char *uplo, char *trans, char *diag, int *n, 
	int *nrhs, double_cx *a, int *lda, double_cx *b, 
	int *ldb, int *info);
 
int ztzrqf_(int *m, int *n, double_cx *a, 
	int *lda, double_cx *tau, int *info);
 
int ztzrzf_(int *m, int *n, double_cx *a, 
	int *lda, double_cx *tau, double_cx *work, int *lwork,
	 int *info);
 
int zung2l_(int *m, int *n, int *k, 
	double_cx *a, int *lda, double_cx *tau, double_cx *
	work, int *info);
 
int zung2r_(int *m, int *n, int *k, 
	double_cx *a, int *lda, double_cx *tau, double_cx *
	work, int *info);
 
int zungbr_(char *vect, int *m, int *n, int *k, 
	double_cx *a, int *lda, double_cx *tau, double_cx *
	work, int *lwork, int *info);
 
int zunghr_(int *n, int *ilo, int *ihi, 
	double_cx *a, int *lda, double_cx *tau, double_cx *
	work, int *lwork, int *info);
 
int zungl2_(int *m, int *n, int *k, 
	double_cx *a, int *lda, double_cx *tau, double_cx *
	work, int *info);
 
int zunglq_(int *m, int *n, int *k, 
	double_cx *a, int *lda, double_cx *tau, double_cx *
	work, int *lwork, int *info);
 
int zungql_(int *m, int *n, int *k, 
	double_cx *a, int *lda, double_cx *tau, double_cx *
	work, int *lwork, int *info);
 
int zungqr_(int *m, int *n, int *k, 
	double_cx *a, int *lda, double_cx *tau, double_cx *
	work, int *lwork, int *info);
 
int zungr2_(int *m, int *n, int *k, 
	double_cx *a, int *lda, double_cx *tau, double_cx *
	work, int *info);
 
int zungrq_(int *m, int *n, int *k, 
	double_cx *a, int *lda, double_cx *tau, double_cx *
	work, int *lwork, int *info);
 
int zungtr_(char *uplo, int *n, double_cx *a, 
	int *lda, double_cx *tau, double_cx *work, int *lwork,
	 int *info);
 
int zunm2l_(char *side, char *trans, int *m, int *n, 
	int *k, double_cx *a, int *lda, double_cx *tau, 
	double_cx *c__, int *ldc, double_cx *work, int *info);
 
int zunm2r_(char *side, char *trans, int *m, int *n, 
	int *k, double_cx *a, int *lda, double_cx *tau, 
	double_cx *c__, int *ldc, double_cx *work, int *info);
 
int zunmbr_(char *vect, char *side, char *trans, int *m, 
	int *n, int *k, double_cx *a, int *lda, double_cx 
	*tau, double_cx *c__, int *ldc, double_cx *work, int *
	lwork, int *info);
 
int zunmhr_(char *side, char *trans, int *m, int *n, 
	int *ilo, int *ihi, double_cx *a, int *lda, 
	double_cx *tau, double_cx *c__, int *ldc, double_cx *
	work, int *lwork, int *info);
 
int zunml2_(char *side, char *trans, int *m, int *n, 
	int *k, double_cx *a, int *lda, double_cx *tau, 
	double_cx *c__, int *ldc, double_cx *work, int *info);
 
int zunmlq_(char *side, char *trans, int *m, int *n, 
	int *k, double_cx *a, int *lda, double_cx *tau, 
	double_cx *c__, int *ldc, double_cx *work, int *lwork,
	 int *info);
 
int zunmql_(char *side, char *trans, int *m, int *n, 
	int *k, double_cx *a, int *lda, double_cx *tau, 
	double_cx *c__, int *ldc, double_cx *work, int *lwork,
	 int *info);
 
int zunmqr_(char *side, char *trans, int *m, int *n, 
	int *k, double_cx *a, int *lda, double_cx *tau, 
	double_cx *c__, int *ldc, double_cx *work, int *lwork,
	 int *info);
 
int zunmr2_(char *side, char *trans, int *m, int *n, 
	int *k, double_cx *a, int *lda, double_cx *tau, 
	double_cx *c__, int *ldc, double_cx *work, int *info);
 
int zunmr3_(char *side, char *trans, int *m, int *n, 
	int *k, int *l, double_cx *a, int *lda, double_cx 
	*tau, double_cx *c__, int *ldc, double_cx *work, int *
	info);
 
int zunmrq_(char *side, char *trans, int *m, int *n, 
	int *k, double_cx *a, int *lda, double_cx *tau, 
	double_cx *c__, int *ldc, double_cx *work, int *lwork,
	 int *info);
 
int zunmrz_(char *side, char *trans, int *m, int *n, 
	int *k, int *l, double_cx *a, int *lda, double_cx 
	*tau, double_cx *c__, int *ldc, double_cx *work, int *
	lwork, int *info);
 
int zunmtr_(char *side, char *uplo, char *trans, int *m, 
	int *n, double_cx *a, int *lda, double_cx *tau, 
	double_cx *c__, int *ldc, double_cx *work, int *lwork,
	 int *info);
 
int zupgtr_(char *uplo, int *n, double_cx *ap, 
	double_cx *tau, double_cx *q, int *ldq, double_cx *
	work, int *info);
 
int zupmtr_(char *side, char *uplo, char *trans, int *m, 
	int *n, double_cx *ap, double_cx *tau, double_cx *c__,
	 int *ldc, double_cx *work, int *info);

#endif /* __CLAPACK_H */
