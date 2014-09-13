#include "mp_mat.h"
#include "lapack.h"

template <>
void dgetrf(mp_mat<double>& A, vector<int> &ipiv) {
     int info = 1;

     if ((int)ipiv.size() != A.n)
	 throw gen_err("dgetrf error: ipiv dimension");
     dgetrf(A.m, A.n, A.p, A.m, &ipiv[0], &info);

     if (info != 0)
 	throw gen_err("info != 0  in dgetrf\n");
}


template <>
void dgetrs(mp_mat<double>& A, vector<int> &ipiv, double *b,
	    int nrhs, char trans)
{
    if (A.m != A.n || A.n != (int)ipiv.size())
	throw gen_err("dgetrs error\n");
    int info = 1;

    dgetrs(trans, A.m, nrhs, A.p, A.m, &ipiv[0], b, A.m, &info);

    if (info != 0)
	throw gen_err("info != 0  in dgetrs\n");
}




// void dgemv(char trans, double alpha, mp_mat<double>& A,
// 	   double *x, double beta, double *y) {
//     int m = A.m;
//     int n = A.n;
//     dgemv(trans, m, n, alpha, A.p, m, x, 1, beta, y, 1);
// }


// void mydgemv(double alpha, mp_mat<double>& A, double *x, double beta,
// 	     double *y) {
//     int m = A.m;
//     int n = A.n;
//     dgemv('N', m, n, alpha, A.p, m, x, 1, beta, y, 1);
// }


// mp_mat<double> inverse(mp_mat<double>& A) {
//     if (A.m != A.n || A.m <= 0)
// 	throw gen_err("mp_mat::inverse() needs square matrix");

//     int n = A.n;
//     mp_mat<double> tmp1;
//     mp_mat<double> tmp2(n,n,0.0);
//     vector<int> ipiv(n);
//     tmp1 = A;

//     for (int i=0; i<n; i++)
// 	tmp2(i,i) = 1.0;

//     dgetrf(tmp1, &ipiv[0]);
//     dgetrs(tmp1, &ipiv[0], tmp2.p, n);
//     return tmp2;
// }

template<> void
dgemm(char trA, char trB, double alpha, mp_mat<double> &A,
      mp_mat<double> &B, double beta, mp_mat<double> &C) {
    int m, n, k, k1;
    if (trA == 'N') {
	m = A.m;
	k = A.n;
    } else {
	m = A.n;
	k = A.m;
    }
    if (trB == 'N') {
	k1 = B.m;
	n = B.n;
    } else {
	k1 = B.n;
	n = B.m;
    }
    if (k1 != k) throw gen_err("dimension error in dgemm");
    dgemm(trA, trB, m, n, k, alpha, &A(0,0), A.m, &B(0,0),
	  B.m, beta, &C(0,0), C.m);
}

template<> void
dgemv(char trans, double a, mp_mat<double> &A, double *x, int incx,
      double b, double *y, int incy)
{
    dgemv(trans, A.m, A.n, a, A.p, A.m, x, incx, b, y, incy);
}

template<>
vector<double> dgesvd(char jobu, char jobv, mp_mat<double>& A,
		      void *UU, void *VV)
{
    // A=USV^T, return diag(S)
    // jobu,v = 'N', 'O', 'S', 'A'  (none, overwrite, compact, all)
    
    mp_mat<double> *U = (mp_mat<double> *) UU;
    mp_mat<double> *VT = (mp_mat<double> *) VV;
    
    if ((jobu == 'S' || jobu == 'A') && UU==NULL ||
	(jobv == 'S' || jobv == 'A') && VV==NULL)
	throw gen_err("dgesvd input error");

    int info=-1;
    int min_nm = (A.m < A.n) ? A.m : A.n;
    vector<double> svals(min_nm,0.0);
    
    double *Up = NULL;
    double *VTp = NULL;
    int ldu = 1;
    int ldv = 1;
    
    if (jobu == 'S' || jobu == 'A') {
	ldu = A.m;
	U->resize( A.m, ((jobu=='A') ? A.m : min_nm) );
	U->reset_values(0.0);
	Up = U->p;
    }
    if (jobv == 'S' || jobv == 'A') {
	ldv = (jobv=='A') ? A.n : min_nm;
	VT->resize( ldv , A.n );
	VT->reset_values(0.0);
	VTp = VT->p;
    }
    int nzr=A.m, nzc=A.n;

    if (jobu != 'A')
	for ( ; nzr>nzc; nzr--) {
	    int j;
	    for (j=0; j<A.n; j++)
		if (A(nzr-1,j) != 0.0) break;
	    if (j<A.n && A(nzr-1,j) != 0.0) break;
	}
    if (jobv != 'A')
	for ( ; nzc>nzr; nzc--) {
	    int i;
	    for (i=0; i<A.m; i++)
		if (A(i,nzc-1) != 0.0) break;
	    if (i<A.m && A(i,nzc-1) != 0.0) break;
	}
    // printf("%d %d %d %d\n", A.m, A.n, nzr, nzc);

    dgesvd(jobu, jobv, nzr, nzc, A.p, A.m,
	   &svals[0], Up, ldu, VTp, ldv, &info);
    
    if (info != 0)
	throw gen_err("dgesvd error");
    
    return svals;
}

template<>
vector<double> dsyev(char jobZ, mp_mat<double>& A)
{
    // A = QSQ^T, return diag(S), overwrite A with Q if jobZ != 'N'
    int info=-1;
    if (A.m != A.n)
	throw gen_err("dsyev input error");
    vector<double> s(A.m,0.0);
    int lda = A.m;
    dsyev(jobZ, 'L', A.m, A.p, lda, &s[0], &info);
    if (info != 0)
	throw gen_err("dsyev error");
    
    return s;
}
