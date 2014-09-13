// ./nlEig 100 50 3 blah > output
//  grep amp output > amp
//  gnuplot
//     plot 'amp' u 2:3 w lp
//     plot 'blah' w lp

#include "levmar.h"


#ifndef STRUCT_CONSTANTS
#define STRUCT_CONSTANTS
struct constants {
    Real pi;
    Real pi2;
    constants() {
	pi = 4.0*atan(Real(1.0));
	pi2 = 2.0*pi;
    }
};
#endif


class nlEig : public levmar, public constants {
public:
    // x = lam,u[1],u[2],..,u[n-1]
    // r = u[n/2]-a, Delta u + lam*sin(u)
    Real amp;
    Real dx;
    nlEig(int n_, Real a);
    void compute_r();
    void compute_J();
};


nlEig::nlEig(int n_, Real a) : levmar(n_,n_), amp(a), dx(pi/n)
{
    if (n%2) throw gen_err("need n even");
    if (a<=0.0) throw gen_err("need amp>0");

    x[0] = 1.0; // small amplitude lambda
    
    vector<Real> &u = x; // alias the name
    for (int i=1; i<n; i++)
	u[i] = amp*sin(i*dx);
}


void nlEig::compute_r()
{
    vector<Real> &u = x;
    Real lam = x[0];
    Real dx2 = dx*dx;

    r[ 0 ] = u[n/2] - amp;
    r[ 1 ] = (  0.0  - 2*u[ 1 ] + u[2] )/dx2 + lam*sin(u[1]);
    r[n-1] = (u[n-2] - 2*u[n-1] + 0.0  )/dx2 + lam*sin(u[n-1]);
    for (int i=2; i<n-1; i++)
	r[i] = ( u[i-1] - 2*u[i] + u[i+1] )/dx2 + lam*sin(u[i]);
}


void nlEig::compute_J()
{
    vector<Real> &u = x;
    Real fac = 1.0/(dx*dx);

    J.reset_values(0.0);
    J(0,n/2) = 1;

    for (int i=1; i<n; i++) {
	J(i,0) = sin(u[i]); // lambda
	if (i>1) J(i,i-1) = fac;
	J(i,i) = -2*fac + x[0]*cos(u[i]);
	if (i<n-1) J(i,i+1) = fac;
    }
}


int main(int argc, char *argv[])
{
    if (argc != 5)
	throw gen_err("usage: ./nlEig n steps amp_max fname");

    int n, steps;
    Real amp_max;

    sscanf(argv[1], "%d", &n);
    sscanf(argv[2], "%d", &steps);
    sscanf(argv[3], "%lf", &amp_max);

    Real amp0 = amp_max / steps;
    nlEig  z(n, amp0);

    z.chkder(5.0e-6);

    for (int i=1; i<=steps; i++) {
	z.amp = amp0*i;
	z.reset_delta();
	z.solve();
	printf("amp_lam  %23.17g  %23.17g\n", z.amp, z.x[0]); // amp,lam
    }
    FILE *fp = fopen(argv[4],"w");
    fprintf(fp, "%23.17g  %23.17g\n", 0.0, 0.0);
    for (int i=1; i<n; i++)
	fprintf(fp, "%23.17g  %23.17g\n", z.pi*i/n, z.x[i]);
    fprintf(fp, "%23.17g  %23.17g\n", z.pi, 0.0);
    fclose(fp);

    return 0;
}
