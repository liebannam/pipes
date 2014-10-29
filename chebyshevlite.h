#ifndef CHEBYSHEVLITE_H
#define CHEBYSHEVLITE_H
#include <math.h>
#include "globals.h"
#include <iostream>
#include "real_def.h"

//A handful of Chebyshev routines (adapted from Jon's chebyshev class)
//
////fill x with n standard chebyshev nodes on [-1,1] 
template<typename T>
void getChebNodes(T *xx, int n)
{
//	Real pi = 4.0*atan(1.);
	for (int i=1; i<n; i++)
	xx[i] = cos( PI*(n-i)/n );
    	xx[0] = -1.0;
    	xx[n] = 1.0;
}

//use 2-term recursion to evaluate Chebyshev polynomial 
//defined over x points in [-1,1] (adapted from Jon's chebyshev class)
//i.e. f(x) = sum_k a_k T_k(x) 
//
template<typename T, typename TT>
TT ChebEval(T *alpha, int n, TT x)
{
	
	double b0, b1, tmp;
	TT u;
	b1 = alpha[n];
	b0 = 2*x*b1 + alpha[n-1];
	for (int i=n-2; i>0; i--) 
	{
		tmp = 2*x*b0 + alpha[i] - b1;
		b1 = b0;
		b0 = tmp;
	}
	u = x*b0 + alpha[0] - b1;
	return u;
}
#endif
