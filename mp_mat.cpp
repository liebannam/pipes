#include "mp_mat.h"
#if !MPFR
#include "fftw3.h"
#endif

template <class T>
mp_mat<T>::mp_mat() :
    m(0), n(0)
{
#if MPFR
    p = new T[m*n]; // so that p can be freed when destroyed
#else
    p = (T *) fftw_malloc(m*n*sizeof(T));
#endif

    num_handles = new int;
    *num_handles = 1;
}

template <class T>
mp_mat<T>::mp_mat(int m_, int n_) : // uninitialized
    m(m_), n(n_)
{
    if (m<0 || n<0)
	throw gen_err("m<0 or n<0 in mp_mat");
#if MPFR
    p = new T[m*n];
#else
    p = (T *) fftw_malloc(m*n*sizeof(T));
#endif

    num_handles = new int;
    *num_handles = 1;
}

#if !MPFR
template <class T>
mp_mat<T>::mp_mat(int m_, int n_, T val) :
    m(m_), n(n_)
{
    if (m<0 || n<0)
	throw gen_err("m<0 or n<0 in mp_mat");
    const int mn = m*n;

    // p = new T[m*n];
    p = (T *) fftw_malloc(m*n*sizeof(T));
    for (int i = 0; i<mn; i++)
	p[i] = val;

    num_handles = new int;
    *num_handles = 1;
}
#endif

template <class T>
mp_mat<T>::mp_mat(int m_, int n_, const char *c) :
    m(m_), n(n_)
{
    istringstream s(c);
    const int mn = m*n;
    string tmp;

    if (m<0 || n<0)
	throw gen_err("error in mp_mat constructor");

#if MPFR
    p = new T[mn];
#else
    p = (T *) fftw_malloc(mn*sizeof(T));
#endif

    int k=0;
    try {
	for (int i=0; i<m; i++)
	    for (int j=0; j<n; j++) {
		// if you run out of data, repeat the last one
		if (s.good()) {
		    k++;
		    s >> tmp;
		}
		istringstream tmp1(tmp);
		tmp1 >> p[i+j*m];
	    }
    }
    catch (...) {
	printf("error in mp_mat(m,n,const char *)");
	exit(1);
    }
    if (k != m*n) {
	if (k==1)
	    printf("warning: only %d object read into %dx%d = %d matrix\n",
		   k, m, n, m*n);
	else
	    printf("warning: only %d objects read into %dx%d = %d matrix\n",
		   k, m, n, m*n);
    }
    num_handles = new int;
    *num_handles = 1;
}


template <class T>
mp_mat<T>::mp_mat(int m_, int n_, istream &s) :
    m(m_), n(n_)
{
    const int mn = m*n;
    string tmp;

    if (m<0 || n<0)
	throw gen_err("error in mp_mat constructor");

#if MPFR
    p = new T[mn];
#else
    p = (T *) fftw_malloc(mn*sizeof(T));
#endif

    int k=0;
    try {
	for (int i=0; i<m; i++)
	    for (int j=0; j<n; j++) {
		// if you run out of data, repeat the last one
		if (s.good()) {
		    s >> tmp;
		    k++;
		}
		istringstream tmp1(tmp);
		tmp1 >> p[i+j*m];
	    }
    }
    catch (...) {
	printf("error in mp_mat(m,n,const char *)");
	exit(1);
    }
    if (k != m*n) {
	if (k==1)
	    printf("warning: only %d object read into %dx%d = %d matrix\n",
		   k, m, n, m*n);
	else
	    printf("warning: only %d objects read into %dx%d = %d matrix\n",
		   k, m, n, m*n);
    }
    num_handles = new int;
    *num_handles = 1;
}

template <class T>
mp_mat<T>::~mp_mat()
{
    if (num_handles) {
	(*num_handles)--;
	if (*num_handles == 0) {
#if MPFR
	    delete[] p;
#else
	    fftw_free(p);
#endif
	    delete num_handles;
	}
    }
}

template <class T>
mp_mat<T>::mp_mat(const mp_mat<T> &A) :
    num_handles(A.num_handles), m(A.m), n(A.n), p(A.p)
{
    if (num_handles) {
	// shallow copy
	(*num_handles)++;
    }
}

template <class T>
mp_mat<T>& mp_mat<T>::copyS(const mp_mat<T> &A)
{
    if (this != &A) {
	if (num_handles) {
	    (*num_handles)--;
	    if (*num_handles == 0) {
#if MPFR
		delete[] p;
#else
		fftw_free(p);
#endif
		delete num_handles;
	    }
	}
	num_handles = A.num_handles;
	m = A.m;
	n = A.n;
	p = A.p;
	if (num_handles)
	    (*num_handles)++;
    }
    return *this;
}

template <class T>
mp_mat<T>& mp_mat<T>::copyD(const mp_mat<T> &A)
{
    if (this != &A) { // beware of self-assignment
	const int mn = A.m*A.n;
	resize( A.m , A.n );
	for (int i = 0; i<mn; i++)
	    p[i] = A.p[i];
    }
    return *this;
}


#if 0  // shallow copy

template <class T>
mp_mat<T>& mp_mat<T>::operator=(const mp_mat<T> &A)
{
    if (this != &A) {
	if (num_handles) {
	    (*num_handles)--;
	    if (*num_handles == 0) {
#if MPFR
		delete[] p;
#else
		fftw_free(p);
#endif
		delete num_handles;
	    }
	}
	num_handles = A.num_handles;
	m = A.m;
	n = A.n;
	p = A.p;
	if (num_handles)
	    (*num_handles)++;
    }
    return *this;
}

#else  // deep copy

template <class T>
mp_mat<T>& mp_mat<T>::operator=(const mp_mat<T> &A)
{
    if (this != &A) { // beware of self-assignment
	const int mn = A.m*A.n;
	resize( A.m , A.n );
	for (int i = 0; i<mn; i++)
	    p[i] = A.p[i];
    }
    return *this;
}

#endif

template <class T>
void mp_mat<T>::resize_and_share_data(int m1, int n1, T *p1)
{
    if (p1 != p) {
	if (num_handles) {
	    (*num_handles)--;
	    if (*num_handles == 0) {
#if MPFR
		delete[] p;
#else
		fftw_free(p);
#endif
		delete num_handles;
	    } 
	}
	num_handles = NULL; // turn off freeing memory
	m = m1;
	n = n1;
	p = p1;
    } else {
	if (m*n != m1*n1)
	    throw gen_err("mp_mat::resize_and_share_data called incorrectly");
	m = m1;
	n = n1;
    }
}


template<class T>
void mp_mat<T>::reset_values(T val)
{
    const int mn = m*n;
    for (int i = 0; i<mn; i++)
	p[i] = val;
}

template<class T>
void mp_mat<T>::dump(const char *fname, int prec) {
    ofstream fp(fname);

    fp.precision(prec);
    mp_mat<T> &A = *this;
    fp << "% load " << fname << "; " << fname <<
	" = reshape(" << fname << "," <<
	m << "," << n << ")" << endl;
    for (int j=0; j<n; j++) {
	for (int i=0; i<m; i++)
	    fp << A(i,j) << endl;
	fp << endl << endl;
    }
}

template<class T>
vector<T> mp_mat<T>::extract_column(int j)
{
    vector<T> col(m);
    if (j<0 || j>=n)
	throw gen_err("mp_mat::extract_column called incorrectly");
    for (int i=0; i<m; i++)
	col[i] = p[i + j*m];
    return col;
}

template<class T>
void mp_mat<T>::resize(int m1, int n1, T val)
{
    resize(m1,n1);
    reset_values(val);
}

template <class T>
void mp_mat<T>::resize(int m1, int n1)
{
    if (m*n == m1*n1) {
	m = m1;
	n = n1;
    } else {
	if (num_handles) {
	    (*num_handles)--;
	    if (*num_handles == 0) {
#if MPFR
		delete[] p;
#else
		fftw_free(p);
#endif
	    } else {
		num_handles = new int;
	    }
	} else {
	    num_handles = new int;
	}	    
	*num_handles = 1;
	m = m1;
	n = n1;
#if MPFR
	p = new T[m*n];
#else
	p = (T *) fftw_malloc(m*n*sizeof(T));
#endif
    }
}


template <class T>
void mp_mat<T>::operator*=(T a)
{
    for (int i=m*n-1; i>=0; i--)
	p[i] *= a;
}

template <class T>
void mp_mat<T>::operator+=(const mp_mat<T> &A) {
    if (m != A.m || n != A.n)
	throw gen_err("dimensions don't match in mp_mat::operator+=");
    int mn = m*n;
    for (int i=0; i<mn; i++)
	p[i] += A.p[i];
}

template <class T>
void mp_mat<T>::operator-=(const mp_mat<T> &A) {
    if (m != A.m || n != A.n)
	throw gen_err("dimensions don't match in mp_mat::operator-=");
    int mn = m*n;
    for (int i=0; i<mn; i++)
	p[i] -= A.p[i];
}

template <class T>
T mp_mat<T>::inner_prod(const vector<T>& u, const vector<T>& v)
{
    if (m != (int) u.size() || n != (int) v.size())
	throw gen_err("mp_mat::inner_product error");
    T sum = T(0.0);
    T tmp, tmp1;
    for (int i=0; i<m; i++) {
	tmp = 0.0;
	for (int j=0; j<n; j++) {
	    tmp1 = p[i+j*m]*v[j]; // avoid temporaries in mpfr_class
	    tmp += tmp1;
	}
	tmp1 = tmp*u[i];
	sum += tmp1;
    }
    return sum;
}


template <class T>
void dump_vec(const vector<T> &x, const char *name) {
    int n = x.size();
    ofstream os(name);
    for (int i=0; i<n; i++)
	os << x[i] << endl;
}

template <class T>
vector<T> operator*(const vector<T> &u, T a) {
    vector<T> ret(u);
    int e = (int) ret.size();
    for (int i=0; i<e; i++)
	ret[i] *= a;
    return ret;
}

template <class T>
T operator*(const vector<T> &u, const vector<T> &v) {
    if (u.size() != v.size())
	throw gen_err("u*v needs vectors to have same size");
    if (u.size() == 0u)
	return T(0.0);
    T ret(0.0);
    typename vector<T>::const_iterator p = u.begin();
    typename vector<T>::const_iterator q = v.begin();
    typename vector<T>::const_iterator e = u.end();
    for (; p!=e; p++, q++)
 	ret += (*p)*(*q);
    return ret;
}

template <class T>
vector<T> operator-(const vector<T> &u) {
    vector<T> v( u.begin(), u.end() );
    typename vector<T>::iterator p = v.begin();
    typename vector<T>::iterator e = v.end();
    for (; p!=e; p++)
 	*p = -(*p);
    return v;
}


template <class T>
vector<T> operator+(const vector<T> &u, const vector<T> &v) {
    if (u.size() != v.size())
	throw gen_err("sizes wrong in vector<T>::operator+");
    vector<T> ret;
    ret.reserve(u.size());
    typename vector<T>::const_iterator p = u.begin();
    typename vector<T>::const_iterator q = v.begin();
    typename vector<T>::const_iterator e = u.end();
    for (; p!=e; p++, q++)
 	ret.push_back(*p + *q);
    return ret;
}

template <class T>
vector<T> operator-(const vector<T> &u, const vector<T> &v) {
    if (u.size() != v.size())
	throw gen_err("sizes wrong in vector<T>::operator-");
    vector<T> ret;
    ret.reserve(u.size());
    typename vector<T>::const_iterator p = u.begin();
    typename vector<T>::const_iterator q = v.begin();
    typename vector<T>::const_iterator e = u.end();
    for (; p!=e; p++, q++)
 	ret.push_back(*p - *q);
    return ret;
}


template <class T>
mp_mat<T> operator*(const mp_mat<T> &A, const mp_mat<T> &B) {
    if (A.n != B.m)
	throw gen_err("matrix dimension mismatch, mp_mat::operator*");

    int m=A.m, K=A.n, n=B.n;
    mp_mat<T> C(m, n);

    for (int j=0; j<n; j++)
	for (int i=0; i<m; i++) {
	    T tmp(0);
	    for (int k=0; k<K; k++)
	        tmp += A(i,k)*B(k,j);
	    C(i,j) = tmp;
	}
    return C;
}


template <class T>
mp_mat<T> operator+(const mp_mat<T> &A, const mp_mat<T> &B) {
    mp_mat<T> C(A.m, A.n);
    if (A.m != B.m  ||  A.n != B.n)
	throw gen_err("dimensions don't match in mp_mat::operator+");
    const int mn = A.m * A.n;
    for (int i=0; i<mn; i++)
	C.p[i] = A.p[i] + B.p[i];
    return C;
}


template <class T>
mp_mat<T> operator-(const mp_mat<T> &A, const mp_mat<T> &B) {
    mp_mat<T> C(A.m, A.n);
    if (A.m != B.m  ||  A.n != B.n)
	throw gen_err("dimensions don't match in mp_mat::operator+");
    const int mn = A.m * A.n;
    for (int i=0; i<mn; i++)
	C.p[i] = A.p[i] - B.p[i];
    return C;
}


template <class T>
mp_mat<T> operator*(const mp_mat<T> &A, T a) {
    mp_mat<T> C; // = A would be shallow copy
    C = A;
    C *= a;
    return C;
}


template <class T>
mp_mat<T> operator*(T a, const mp_mat<T> &A) {
    return A*a; // copy constructor is cheap
}


template <class T>
mp_mat<T> trans(const mp_mat<T> &A) {
    mp_mat<T> B(A.n, A.m);
    transpose(A, B);
    return B;
}


template <class T>
vector<T> operator*(const mp_mat<T> &A, const vector<T> &x)
{
    if ( A.n != (int)x.size() )
 	throw gen_err("dimensions of A and x don't match in (mp_mat) A*x");

    int m=A.m, n=A.n;
    vector<T> b;
    b.reserve(m);

    for (int i=0; i<m; i++) {
	T tmp(0);
	for (int j=0; j<n; j++)
	    tmp += A(i,j)*x[j];
	b.push_back(tmp);
    }
    return b;
}


template <class T>
vector<T> operator*(const vector<T> &x, const mp_mat<T> &A) {
    if ( A.m != (int)x.size() )
 	throw gen_err("dimensions of A and x don't match in (mp_mat) A*x");
    int m=A.m, n=A.n;
    vector<T> b;
    b.reserve(n);

    for (int j=0; j<n; j++) {
	T tmp(0);
	for (int i=0; i<m; i++)
	    tmp += x[i]*A(i,j);
	b.push_back(tmp);
    }
    return b;
}

template <class T>
void transpose(const mp_mat<T>& A, mp_mat<T>& B)
{
    if (A.p == B.p)
	throw gen_err("error: transpose of matrices sharing data");
	
    if (B.m != A.n  ||  B.n != A.m )
	B.resize(A.n, A.m);

    for (int i=0; i<A.m; i++)
	for (int j=0; j<A.n; j++)
	    B(j,i) = A(i,j);
}

template<class T>
T compute_norm(vector<T> &x) {
    T ret = T(0.0);
    int n = x.size();
    for (int i=0; i<n; i++)
	ret += x[i]*x[i];
    return sqrt(ret);
}

template<class T>
void scale_vec(vector<T> &x, const vector<T> &D) {
    int n = (int) x.size();
    if (n != (int)D.size())
	throw gen_err("error in scale_vec");
    for (int i=0; i<n; i++)
	x[i] *= D[i];
}

template<class T>
void scale_mat(mp_mat<T> &A, const vector<T> &s)
{
    int m = A.m;
    int n = A.n;
    if (n != (int) s.size())
	throw gen_err("error in scale_mat");
    for (int j=0; j<n; j++)
	for (int i=0; i<m; i++)
	    A(i,j) *= s[j];
}

template<class T>
void scale_mat(const vector<T> &s, mp_mat<T> &A)
{
    int m=A.m;
    int n=A.n;
    if (m != (int) s.size())
	throw gen_err("error in scale_mat");
    for (int j=0; j<n; j++)
	for (int i=0; i<m; i++)
	    A(i,j) *= s[i];
}

template<class T>
void sandwichL(char trans, mp_mat<T>& X, mp_mat<T>& K, mp_mat<T>& W)
{
    if (K.m != K.n ||
	(trans == 'N' && K.n != W.m) ||  // X = WT*K*W
	(trans != 'N' && W.n != K.n))    // X = W*K*WT
	throw gen_err("error in mp_mat::sandwichL");

    int n = (trans == 'N') ? W.n : W.m;
    int k = K.m;
    X.resize(n,n);
    vector<T> tmp(k);
    T sum;

    if (trans == 'N') {
	// X = WT*K*W
	for (int j=0; j<n; j++) {
	    for (int a=0; a<k; a++) {
		sum=0.0;
		for (int b=0; b<a; b++)
		    sum += K(a,b)*W(b,j);
		for (int b=a; b<k; b++)
		    sum += K(b,a)*W(b,j);
		tmp[a] = sum;
	    }
	    for (int i=0; i<j; i++)
		X(i,j) = X(j,i);
	    for (int i=j; i<n; i++) {
		sum=0.0;
		for (int a=0; a<k; a++)
		    sum += tmp[a]*W(a,i);
		X(i,j) = sum;
	    }
	}
    } else {
	// X = W*K*WT
	for (int j=0; j<n; j++) {
	    for (int a=0; a<k; a++) {
		sum=0.0;
		for (int b=0; b<a; b++)
		    sum += K(a,b)*W(j,b);
		for (int b=a; b<k; b++)
		    sum += K(b,a)*W(j,b);
		tmp[a] = sum;
	    }
	    for (int i=0; i<j; i++)
		X(i,j) = X(j,i);
	    for (int i=j; i<n; i++) {
		sum=0.0;
		for (int a=0; a<k; a++)
		    sum += tmp[a]*W(i,a);
		X(i,j) = sum;
	    }
	}
    }
}

#define T double
#include "mp_mat_aux.cpp"
#undef T
#define T complex<double>
#include "mp_mat_aux.cpp"
#undef T

#if DDDD
#include <qd/dd_real.h>
#define T dd_real
#include "mp_mat_aux.cpp"
#undef T
#define T complex<dd_real>
#include "mp_mat_aux.cpp"
#undef T
#endif

#if MPFR
#include "gmpfrxx.h"
#define T mpfr_class
#include "mp_mat_aux.cpp"
#undef T
// for complex<mpfr_class> to work, you have to modify <complex> to avoid ?:
// #define T complex<mpfr_class>
// #include "mp_mat_aux.cpp"
// #undef T
#endif
