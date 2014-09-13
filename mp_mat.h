#ifndef __MP_MAT_H
#define __MP_MAT_H

#include "top.h"
using namespace std;

template<class T> class mp_mat {
private:
    int *num_handles;
    // static vector<mp_mat*> stack_tmp; // assign to them first
public:
    int m, n;
    T *p;

    mp_mat();
    mp_mat(int m, int n);        // uninitialized
    mp_mat(int m, int n, T val); // initialized with val
    mp_mat(int m, int n, const char *s);
    mp_mat(int m, int n, istream &s); // ifstream or istringstream
    ~mp_mat();

    // copy constructor only copies pointer; use vector<mp_mat> carefully
    mp_mat(const mp_mat &);
    mp_mat& operator=(const mp_mat &A); // assignment does deep copy
    
    mp_mat& copyS(const mp_mat &A); // shallow
    mp_mat& copyD(const mp_mat &A); // deep
    void reset_values(T val);
    void resize(int m1, int n1); // uninitialized (m==m1,n==n1: nothing done)
    void resize(int m1, int n1, T val);
    void resize_and_share_data(int m1, int n1, T *p1);
    void operator*=(T a);
    void operator-=(const mp_mat &A);
    void operator+=(const mp_mat &A);

    // A(i,j) = p[i + j*m],  0 <= i < m,  0 <= j < n
    T& operator()(int i, int j) {return p[i+j*m];}
    const T& operator()(int i, int j) const {return p[i+j*m];}
    
    // void dump();
    void dump(const char *fname, int prec=17);

    std::vector<T> extract_column(int j);

    int cols() const { return n; }
    int rows() const { return m; }
    void swapCols(int j, int k) {
	for (int i=0; i<m; i++) {
	    T tmp = p[i+j*m];
	    p[i+j*m] = p[i+k*m];
	    p[i+k*m] = tmp;
	}
    }
    T inner_prod(const vector<T>& u, const vector<T>& v);
}; // end template mp_mat


// void canonicalize(mp_mat<mpq_class> &A);

template <class T> void dump_vec(const vector<T> &x, const char *name);
template <class T> T operator*(const vector<T> &u, const vector<T> &v);
template <class T> vector<T> operator*(const vector<T> &u, T a);
template <class T> vector<T> operator-(const vector<T> &u);
template <class T> vector<T> operator+(const vector<T> &u, const vector<T> &v);
template <class T> vector<T> operator-(const vector<T> &u, const vector<T> &v);
template <class T> mp_mat<T> operator*(const mp_mat<T> &A, const mp_mat<T> &B);
template <class T> mp_mat<T> operator+(const mp_mat<T> &A, const mp_mat<T> &B);
template <class T> mp_mat<T> operator-(const mp_mat<T> &A, const mp_mat<T> &B);
template <class T> mp_mat<T> operator*(const mp_mat<T> &A, T a);
template <class T> mp_mat<T> operator*(T a, const mp_mat<T> &A);
template <class T> mp_mat<T> trans(const mp_mat<T> &A);
template <class T> vector<T> operator*(const mp_mat<T> &A, const vector<T> &x);
template <class T> vector<T> operator*(const vector<T> &x, const mp_mat<T> &A);
template <class T> void transpose(const mp_mat<T>& A, mp_mat<T>& B);

template<class T> T compute_norm(vector<T> &x);
template<class T> void scale_vec(vector<T> &x, const vector<T> &s);
template<class T> void scale_mat(mp_mat<T> &A, const vector<T> &s);
template<class T> void scale_mat(const vector<T> &s, mp_mat<T> &A);


// jobu,v = 'N', 'O', 'S', 'A'  (none, overwrite, compact, all)
template <class T>
std::vector<T> dgesvd(char jobu, char jobv,
		      mp_mat<T>& A,
		      void *U = NULL,
		      void *VT = NULL);

// jobZ = 'N', 'V'  (compute evecs),  only looks at lower triangle of A
template <class T>
std::vector<T> dsyev(char jobZ, mp_mat<T>& A);

// 'N':  X = WT*K*W,     'T':  X = W*K*WT
// only access/store lower triangle of X and K (assume K symmetric)
template <class T>
void sandwichL(char trans, mp_mat<T>& X, mp_mat<T>& K, mp_mat<T>& W);

template<class T>
void dgemm(char trA, char trB, T a, mp_mat<T> &A, mp_mat<T> &B,
	   T b, mp_mat<T> &C);

template<class T>
void dgemv(char trans, T a, mp_mat<T> &A, T *x, int incx,
	   T b, T *y, int incy);
// DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

template<class T>
void dgetrf(mp_mat<T> &A, vector<int> &ipiv);
template<class T>
void dgetrs(mp_mat<T>& A, vector<int> &ipiv, T *b, int nrhs=1,
	    char trans='N');

#ifndef GEN_ERR
#define GEN_ERR
struct gen_err {
    const char *msg;
    gen_err(const char *p) : msg(p) {
	std::cerr << p << "\n";
    }
};
#endif


#endif // __mp_mat_h
