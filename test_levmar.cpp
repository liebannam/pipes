// ./test_levmar 50 25 | grep delta

#include "levmar.h"

class example : public levmar {
public:
    mp_mat<Real> A;
    vector<Real> r0;
    vector<Real> x0;
    example(int m, int n) :
	levmar(m,n), A(m,n), r0(m), x0(n) {
	for (int j=0; j<n; j++) {
	    for (int i=0; i<m; i++)
		A(i,j) = drand48() - .5;
	    x0[j] = drand48()-.5;
	}
	r0 = A*x0;
    }
    void compute_r() {
	r = A*x - r0;
    }
    void compute_J() {
	copy(A.p, A.p+m*n, J.p);
    }
};

int main(int argc, char *argv[])
{
    if (argc != 3)
	throw gen_err("usage: ./test_levmar m n");

    int m = 100, n = 50;
    sscanf(argv[1], "%d", &m);
    sscanf(argv[2], "%d", &n);

    example a(m,n);
   for(int i =0; i<m; i++){for (int j=0; j<n;j++){printf("A0[%d,%d] = %f\n", i,j,a.A(i,j));} }
   a.set_delta(.6);
   a.compute_r();
   a.compute_J();   
    a.solve();
  //dgemv('T', m, n, Real(1.0), a.J.p, m, &a.r[0], 1, Real(0.0), &a.g[0], 1);
  //  a.compute_f();
   a.compute_g();
    a.dump();

    return 0;
}
