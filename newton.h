

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

template<class F, class dF>
double Newton(F fcn, dF dfcn, double x0, double tol, int maxiter){
	double x = x0;
	int i;
	for(i = 0; i<maxiter; i++){
		if (fabs(dfcn(x0))<1e-16){
			printf("df is singular at x0, aborting Newton's\n");
			exit(1);
		}
		x = x0-fcn(x0)/dfcn(x0);
		if(fabs(x-x0)<tol)
			break;
		x0 = x;
	}
	//if (i == maxiter-1)
	//	printf("Maximum iterations reached\n");
	//else
	//	printf("Converged in %d iterations; f(x) = %f\n", i, fcn(x));
	return x;
	
}
