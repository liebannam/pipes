#ifndef TOP_H

#define TOP_H 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <string>

using namespace std;

// better to use the stl than risk macro clashes
// #define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
// #define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
// #define MAX(X,Y)  ((X) > (Y) ? (X) : (Y))
// #define MIN(X,Y)  ((X) < (Y) ? (X) : (Y))
// #define SQR(X)    ((X)*(X))

//#ifndef PI
//const double    PI = 3.14159265358979323846;
//#endif
const double   PI2 = 6.28318530717958647693;
const double SQRT2 = 1.41421356237309504880;
const double  LOG2 = 0.693147180559945309417;
const double  LOG10 = 2.30258509299404568402;
const double  LOG2PI = 1.83787706640934548356;

#ifndef GEN_ERR
#define GEN_ERR
struct gen_err {
    const char *msg;
    gen_err(const char *p) : msg(p) {
	cout << p << "\n";
    }
};
#endif

typedef std::complex<double> Complex;

class cx_out { // output  z=x+iy  as  x  when  y==0
public:
    double re, im;
    cx_out(Complex z) : re(real(z)), im(imag(z)) { }
};
ostream& operator << (ostream& s, const cx_out& z);

#endif
