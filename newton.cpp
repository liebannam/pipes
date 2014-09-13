


#include <stdio.h>
#include <math.h>

template<class F>
double = 100;
double Newton(){F fcn, F dfcn, double x0, double tol){
	double x = x0;
	int i=0;
	
	for(i = 0; i<maxiter; i++){
		if (dfcn(x0)==0){
			printf("df is singular at x0, aborting Newton's\n");
			exit (1)
		}
		x = x0-fcn(x0)/dfcn(x0);
		if(abs(x-x0)<tol)
			break;
		x0 = x;
	}
	if (i == maxiter-1)
		printf("Maximum iterations reached\n");
	else
		printf("Converged in %d iterations; f(x) = %f\n", i, fcn(x));
	return x;

}
