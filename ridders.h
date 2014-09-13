//#include <cmath>
//#include <cstdio>
//#include <climits>
//#include <cfloat>
#include "top.h"

// modeled on numerical recipes code
template<class F>
double ridders(F fcn, double xl, double xr, int *count,
	     double tol1, double tol2) {
    // notes: error1, error2 never happen (barring atrocious arithmetic or NaN)
    // max iters never reached if 2^-jmax < tol1 (before scaling below)

    int j, TEST = 0;
    const int jmax = 250;
    const int standard = 0;
    double s, xm, xnew=0.0, xold=0.0, xtrial, fl, fr, fm, fnew=0.0, ftrial;
    bool xold_valid = false;
    bool xnew_is_xm = false;
    
    //tol1 *= (xr-xl);
    //tol2 *= (xr-xl);
    fl = fcn(xl);
    fr = fcn(xr);
    if (count) *count = 2;
    
    if (fl==0) return xl;
    if (fr==0) return xr;
    if (fl>0 && fr>0 || fl<0 && fr<0) {
	printf("xl = %f, xr = %f, fl = %f, fr = %f\n",
	       xl, xr, fl, fr);
	throw gen_err("root not bracketed in ridders");
    }

    for (j=1; j<=jmax; j++) {
	
	// get xm
	if (xnew_is_xm) {
	    xm = xnew; fm = fnew;
	} else {
	    xm = .5*(xl+xr); fm = fcn(xm); if (count) ++*count;
	}
	if (fm==0) return xm;
	xnew_is_xm = false;
	
	// get xnew
	s = sqrt(fm*fm-fl*fr);
	if (s<abs(fm)) s = abs(fm); // underflow or non-IEEE roundoff error
	if (fl<fr)
	    xnew = xm - (xm-xl)*fm/s;
	else
	    xnew = xm + (xm-xl)*fm/s;
	if (xnew<xl) xnew=xl; // unusual roundoff error
	if (xnew>xr) xnew=xr;
	
	// print diagnostics if requested
	if (TEST) printf("xl = %f, xnew = %f, xr = %f, xr-xl = %f\n",
			 xl, xnew, xr, xr-xl);
	
	// test if we are done
	if (xold_valid  &&  abs(xnew-xold) <= tol1  &&  xr-xl <= tol2)
	    return xnew;
	xold = xnew; xold_valid = true;
	fnew = fcn(xnew); if (count) ++*count;
	if (fnew == 0.0) return xnew;
	
	// discard half of the interval
	if (fm>0 && fl<0 || fm<0 && fl>0) {xr = xm; fr = fm;}
	if (fm>0 && fr<0 || fm<0 && fr>0) {xl=xm; fl = fm;}
	if (xnew<xl || xnew>xr) {printf("error1 in ridders\n"); exit(1);}
	
	// choose the next interval
	if (fnew>0 && fl<0 || fnew<0 && fl>0) {
	    if (standard || xnew-xl <= 2*(xr-xnew) || xr-xnew==0) {
		xr=xnew; fr=fnew; // usual ridders' method
	    }
	    else {
		xtrial = xr-2*(xr-xnew);
		ftrial = fcn(xtrial); if (count) ++*count;
		if (ftrial==0) return xtrial;
		if (ftrial>0 && fr<0 || ftrial<0 && fr>0) {
		    xl=xtrial; fl=ftrial; // my approach
		    xnew_is_xm = true;
		}
		else {
		    if (TEST) printf("wasted function call\n");
		    xr=xnew; fr=fnew; // usual ridders' method
		}
	    }
	}
	else if (fnew>0 && fr<0 || fnew<0 && fr>0) {
	    if (standard || xr-xnew <= 2*(xnew-xl) || xnew-xl==0) {
		xl=xnew; fl=fnew; // usual ridders' method
	    }
	    else {
		xtrial = xl+2*(xnew-xl);
		ftrial = fcn(xtrial); if (count) ++*count;
		if (ftrial==0) return xtrial;
		if (ftrial>0 && fl<0 || ftrial<0 && fl>0) {
		    xr=xtrial; fr=ftrial; // my approach
		    xnew_is_xm = true;
		}
		else {
		    if (TEST) printf("wasted function call\n");
		    xl=xnew; fl=fnew; // usual ridders' method
		}
	    }
	}
	else {
	    throw gen_err("error2 in ridders; probably fcn returned NaN");
	}
    } // end for j
    printf("max iterations reached in ridders\n");
    return xnew;
} // end ridders


template<class Real, class F>
Real brent(F fcn, Real x1, Real x2, Real fa, Real fb,
	   int *count, Real eps = DBL_EPSILON, Real tol = DBL_EPSILON) {
    int j;
    int b_last = 0;
    const int jmax = 60;
    Real a=x1, b=x2, c=x2, d=0.0, e=0.0, min1, min2;
    // Real fa=fcn(a), fb=fcn(b);
    Real fc=fb, p, q, r, s, tol1, xm;
    
    //tol *= x2-x1;
    if (count) *count = 2;
    
    if (fa>0 && fb>0 || fa<0 && fb<0) {
	printf("a = %s, b= %s, fa = %s, fb = %s\n",
	       str(a), str(b), str(fa), str(fb));
	throw gen_err("root must be bracketed in brent\n");
    }
    for (j=1; j<=jmax; j++) {
	if ((fb>0 && fc>0) || (fb<0 && fc<0)) {
	    c=a; fc=fa; e=d=b-a;
	}
	if (abs(fc)<abs(fb)) {
	    a=b; b=c; c=a; fa=fb; fb=fc; fc=fa; b_last=0;
	}
	if (0) printf("b %17s, c %17s, a %17s, "
		      "b-c %7s\n", str(b,15), str(c,15), str(a,15),
		      str(b-c,5));
	tol1=2.0*eps*abs(b); // +.5*tol;
	xm=.5*(c-b);
	if (abs(xm) <= tol1 || fabs(fb)<=tol) {
	    if (!b_last) {
		fcn(b); // in case fcn causes side effects
		if (count) ++*count;
		// printf("again\n");
	    }
	    return b;
	}
	if (abs(e) >= tol1 && abs(fa)>abs(fb)) {
	    s=fb/fa;
	    if (a==c) {
		p=2.0*xm*s;
		q=1.0-s;
	    } else {
		q=fa/fc;
		r=fb/fc;
		p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
		q=(q-1.0)*(r-1.0)*(s-1.0);
	    }
	    if (p>0.0) q = -q;
	    p=abs(p);
	    min1=3.0*xm*q-abs(tol1*q);
	    min2=abs(e*q);
	    if (2.0*p < (min1<min2 ? min1 : min2)) {
		e=d; d=p/q;
	    } else {
		d=xm; e=d;
	    }
	} else {
	    d=xm; e=d;
	}
	a=b;
	fa=fb;
	if (abs(d)>tol1)
	    b+=d;
	else if (xm >= 0.0)
	    b += abs(tol1);
	else
	    b -= abs(tol1);
	fb=fcn(b);
	b_last = 1;
	if (count) ++*count;
    }
    printf("error: brent's method failed to converge\n");
    return b;
}
