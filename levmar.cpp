#include "levmar.h"

#define DEBUG 0
#define COMPARE_JQ 0
#define CHECK_ORTHOG 0

// note: indices in most comments are 1-based (as in matlab/fortran)


levmar::levmar(int m_, int n_) :
    m(m_), n(n_), 
    x(n,0.0), p(n,0.0), r(m,0.0), g(n,0.0),
    old_x(n,0.0), old_r(m,0.0),
    D(n,1.0), p2(n,0.0), g2(n,0.0), tmp(n,0.0),
#if SCA
    J(0,0), UU(0,0), VT(0,0),
#else
    J(m,n), UU(m,n), VT(n,n),
#endif
    Sig(n,0.0), Lam(n,0.0)
{
    if (n<=0 || m<n) throw gen_err("need m >= n > 0 in levmar");
    set_init_params();
#if SCA
    scaJq = scaUU = scaVT = scaN1 = scaN2 = scaM1 = NULL;
#endif
}


void levmar::set_init_params()
{
    f = 0.0;
    df = 0.0;
    
    norm_g0 = 0.0;
    norm_g = 0.0;
    norm_p2 = 0.0;

    old_f = 0.0;
    old_norm_g = 0.0;
    
    lambda = 0.0;
    verbose = 2;
    update_delta0 = true;

    delta = 1.0;
    delta0 = 1.0;
    delta_max = 1.0e6;
    delta_trigger = -1.0;
    delta_next = delta; // not used unless delta_trigger >= 0.0
    gammaD = 0.0;

    nfev = 0;
    nJev = 0;
    maxfev = 20000;

    num_steps = 0;
    stepsJ = 0;
    max_stepsJ = 10;
    report_r_stride = 10;
    newton_iters = 0;
    
    output_intermediate_successes = 0;

#if DDDD
    f_tol = 1.0e-50; // 1.0e-34
    g_tol = 1.0e-28; // 1.0e-19
    df_tol = 1.0e-5;
    newton_tol = 1.0e-24;
    diag_tol = 1.0e-16;
    sig_tol = 1.0e-26;
#else
    f_tol = 1.0e-25;
    g_tol = 1.0e-15;
    df_tol = 1.0e-5;
    newton_tol = 1.0e-12;
    diag_tol = 1.0e-8;
    sig_tol = 1.0e-13;
#endif

}

#if SCA
void levmar::set_dist_mat_ptrs(dist_mat *a, dist_mat *b, dist_mat *c,
			       dist_mat *d, dist_mat *e, dist_mat *f)
{
	scaJq = a;
	scaUU = b;
	scaVT = c;
	scaN1 = d;
	scaN2 = e;
	scaM1 = f;
	if ( scaJq->m != m  ||  scaJq->n != n ) {
	    char buf[64];
	    sprintf(buf, "dimension error in scalapack pointers: %d  %d,  %d  %d",
		   scaJq->m, m, scaJq->n, n);
	    throw gen_err(buf);
	}
}
#endif

void levmar::set_x(Real *xx)
{
    x.assign(xx,xx+n);
}


void levmar::update_x()
{
    for (int i=0; i<n; i++)
	x[i] += p[i];
}


void levmar::set_gammaD(Real gD)
{
    if ( gD<0.0 || gD>=1 )
	throw gen_err("need 0.0 <= gammaD < 1.0");
    gammaD = gD;
}


void levmar::compute_diag(int nD)
{
    Real minD = 0.0;
    Real fac = 1.0/diag_tol;

    for (int i=0; i<nD; i++) {
	D[i] = 1.0;
	if ( Sig[i] != 0.0  &&  gammaD != 0.0 )
	    D[i] = pow( Sig[i], -gammaD );
	if ( i==0  ||  minD > D[i] )
	    minD = D[i];
    }
    int cnt = 0;
    for (int i=0; i<nD; i++) {
	if (D[i] > fac*minD) {
	    D[i] = fac*minD; // don't let D get too large
	    cnt++;
	}
    }
    if ( cnt>0  &&  verbose>1 )
	printf("diag_tol activated %d times\n", cnt);
}
    

void levmar::dump(FILE *fp) // fixme: needs work. output r?
{
    fprintf(fp,"nfev = %d, nJev = %d\n", nfev, nJev);
    fprintf(fp,"f = %f\n",f);

    compute_g(); // from r and J

    fprintf(fp,"x,g =\n");
    for (int i=0; i<n; i++)
	fprintf(fp,"%4d %23s  %23s\n", i, str(x[i],0), str(g[i],0));
   // J.dump("JJ");
}


void levmar::chkder(Real dx)
{
    // x is the base point
    Real fp;
    Real fq;
    vector<Real> xp(n);
    vector<Real> xq(n);
    vector<Real> rp(m);
    vector<Real> rq(m);

    vector<Real> e(n);
    mp_mat<Real> E(m, n);

    for (int j=0; j<n; j++) {
	for (int i=0; i<n; i++) {
	    int d = (i==j) ? 1 : 0;
	    xp[i] = x[i] + dx*d;
	    xq[i] = x[i] - dx*d;
	}
	x.swap(xp);
	r.swap(rp);
	compute_f();
	fp = f;
	x.swap(xp);
	r.swap(rp);

	x.swap(xq);
	r.swap(rq);
	compute_f();
	fq = f;
	x.swap(xq);
	r.swap(rq);

	e[j] = (fp - fq)/(2.0*dx);
	for (int i=0; i<m; i++)
	    E(i,j) = (rp[i] - rq[i])/(2.0*dx);
    }
    compute_f();
    compute_J();
    compute_g();

    Real max_g = -1.0;
    Real max_gdiff = -1.0;
    for (int i=0; i<n; i++) {
	if (max_g < abs(g[i]))
	    max_g = abs(g[i]);
	if (max_gdiff < abs(g[i] - e[i]))
	    max_gdiff = abs(g[i] - e[i]);
    }
    printf("grad,  (num_grad - grad)\n");
    for (int i=0; i<n; i++)
	printf("%12s  %12s\n", str(g[i]), str(e[i] - g[i]));

    Real max_J = -1;
    Real max_Jdiff = -1;
#if SCA
    int flags[4] = { 0 };
    flags[2] = n;
    flags[1] = 5;
    MPI_Bcast(flags, 4, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    scaJq->consolidate( Jq.p, 0, n );
#endif
    for (int j=0; j<n; j++) {
	for (int i=0; i<m; i++) {
	    if (max_J < abs(J(i,j)))
		max_J = abs(J(i,j));
	    if (max_Jdiff < abs(E(i,j) - J(i,j)))
		max_Jdiff = abs(E(i,j) - J(i,j));
	}
    }
    if (1) {
//	int mc = min(m,20);
	int mc = min(m,101);
	int nc = min(n,10);
	int off1 = 0; // m-mc; // 0 or m-mc
	int off2 = 0; // n-nc; // 0 or n-nc
	printf("Jacobian\n     ");
	for (int i=0; i<nc; i++)
	    printf(" %12d", off2+i);
	printf("\n");
	for (int i=0; i<mc; i++) {
	    printf("%4d ", off1+i);
	    for (int j=0; j<nc; j++)
		printf(" %12s", str(J(off1+i,off2+j)));
	    printf("\n");
	}
	if (0) {
	    printf("numerical Jacobian\n");
	    for (int i=0; i<mc; i++) {
		printf("%4d ", off1+i);
		for (int j=0; j<nc; j++)
		    printf(" %12s", str(E(off1+i,off2+j)));
		printf("\n");
	    }
	} else {
	    printf("numerical Jacobian - Jacobian\n");
	    for (int i=0; i<mc; i++) {
		printf("%4d ", off1+i);
		for (int j=0; j<nc; j++)
		    printf(" %12s", str(E(off1+i,off2+j) - J(off1+i,off2+j)));
		printf("\n");
	    }
	}
    }
    printf("chkder results:\n");
    printf("max abs entry of g: %s\n", str(max_g));
    printf("         g_num - g: %s\n", str(max_gdiff));

    printf("max abs entry of J: %s\n", str(max_J));
    printf("         J_num - J: %s\n", str(max_Jdiff));
}


void levmar::compute_f()
{
    compute_r();
    f = 0.0;
  //  for (int i=0; i<n; i++) I CHANGED THIS MAYBE I'M WRONG!?!?!
    for (int i=0; i<m; i++)
	f += r[i]*r[i];
    f *= 0.5;
}

void levmar::compute_g()
{
#if !SCA
 for(int i =0; i<n*m; i++)
// {cout<<"i = "<<i<< "  and "<<J.p[i]<<endl;} 
//  printf("m = %d, and n = %d\n", m,n);
//  printf("length of r is %d, length of g is %d , size of J is %d \n", sizeof(r)/sizeof(r[0]), sizeof(g)/sizeof(g[0]),sizeof(J.p)/sizeof(J.p[0]));
// dgemv('T', m, n, Real(1.0), J.p, m, &r[0], 1, Real(0.0), &g[0], 1);
  dgemv('T', m, n, Real(1.0), J.p, m, &r[0], 1, Real(0.0), &g[0], 1);
#else
//    throw gen_err("compute_g not implemented for scalapack yet");
#endif
    norm_g = compute_norm(g);
}


void levmar::solve(int skipJ)
{
    nfev = 0;
    nJev = 0;

    num_steps = 0;
    // 2 means set delta0 to delta on 2nd call of check_reduction()
    update_delta0_counter = 2;
    cout<<"Damn."<<endl;
    compute_f(); nfev++;
    old_x = x;

    int flag = 1;
    while (flag) {
	if (!skipJ) {
	    compute_J();
	    nJev++;
	    cout<<nJev<<endl;
	}
	skipJ = 0;
	
    	cout<<"Damn 1."<<endl;
	flag = take_step(); // update x, compute new r
	norm_r = compute_norm(r);
    }
    
    if (report_r_stride>0)
	report_r();
    

}


void levmar::report_r(FILE *fp)
{
    
    fprintf(fp, "\nx =\n");
    for (int i=0; i<n; i++) {
	if (i%8 != 7 && i<n-1)
	    fprintf(fp, "%5d %12s,", i, str(x[i]));
	else
	    fprintf(fp, "%5d %12s\n", i, str(x[i]));
    }
    fflush(fp);
    fprintf(fp, "\nr =\n");
    for (int i=0; i<m; i++) {
	if (i%8 != 7 && i<m-1)
	    fprintf(fp, "%5d %12s,", i, str(r[i]));
	else
	    fprintf(fp, "%5d %12s\n", i, str(r[i]));
    }
    fflush(fp);
    cout<<"fuck?1"<<endl;
}


void reverse(int n, Real *q) {
    for (int i=0; i<n/2; i++) {
	Real tmp = q[i];
	q[i] = q[n-1-i];
	q[n-1-i] = tmp;
    }
}


// make sure f,r,Jq already evaluated at x
int levmar::take_step()
{
    if ( report_r_stride > 0  &&  (num_steps++) % report_r_stride == 0 )
	report_r();

    stepsJ = 0;
    delta_trigger = -1.0;

#if SCA
    int flags[4] = { 0 };
    flags[1] = 2;
    MPI_Bcast(flags, 4, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    scaJq->svd( scaUU, scaVT, &Sig[0] ); // destroys Jq
#else

    UU.resize(m,n);
    VT.resize(n,n);
    Lam.resize(n);
    D.resize(n);
    p2.resize(n);
    copy(J.p, J.p + m*n, UU.p); // Jq remains valid in this version
    if (m<n) throw gen_err("need m>=n in levmar::take_step");
    cout<<"Damn! need m>=n, we have m = "<<m<<"n = "<<n<<endl;
    Sig = dgesvd('O', 'S', UU, NULL, &VT);  //THIS IS THE PROBLEM LINE --throws segfault here...
#endif

    for (int i=0; i<n; i++)
	if ( Sig[i] < sig_tol*Sig[0] )
	    Sig[i] = 0.0;

    compute_diag(n);
    
    cout<<"Damn 1.8"<<endl;
    if (verbose>1) {
	printf("\n%12s %12s\n", "Sig", "D");
	for (int j=0; j<n; j++)
	    printf("%10d %12s %12s\n", j, str(Sig[j]), str(D[j]));
	printf("%12s %12s\n\n", "Sig", "D");
    }
    
    for (int i=0; i<n; i++)
	Lam[n-1-i] = D[i]*Sig[i] * D[i]*Sig[i];
    
    while (nfev < maxfev) {
#if SCA
	flags[1] = 3;
	MPI_Bcast(flags, 4, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	scaM1->distribute(&r[0]);
	scaUU->apply_to_vec( *scaM1, *scaN1, 'T', 1.0, 0.0 );
	scaN1->consolidate(&g2[0]);
#else
	// g2 = r*UU;
	dgemv('T', m, n, Real(1.0), UU.p, m, &r[0], 1, Real(0.0), &g2[0], 1);	
#endif
	scale_vec(g2, Sig);
	norm_g = compute_norm(g2); // |g| computed before g2 rescaled by D
	scale_vec(g2, D);
	reverse(n,&g2[0]);

	compute_p2(); // also df, norm_p2
	reverse(n,&p2[0]);
	p = p2; // p2 and D have size n
	scale_vec(p, D);
#if SCA
	scaN1->distribute(&p[0]);
	scaVT->apply_to_vec( *scaN1, *scaN2, 'T', 1.0, 0.0 );
	scaN2->consolidate(&p[0]);
#else
	// p = p*VT;
	tmp = p;
	dgemv('T', n, n, Real(1.0), VT.p, n, &tmp[0], 1, Real(0.0), &p[0], 1);
#endif

    	cout<<"Damn 2."<<endl;
	if (stepsJ==0) {
	    norm_g0 = norm_g; // norm_g0 part of roundoff regime decision
	    if (verbose)
		printf("norm_g0 = %s\n", str(norm_g0));
	} else {
	    if ( f >= f_tol && delta >= g_tol && norm_g0 >= g_tol
		 && abs(df) < df_tol*f ) {
		// not in roundoff regime, stepsJ>0, df is tiny. time to recompute J
		if (verbose>1)
		    printf("short circuiting stepsJ. |df| < df_tol*f\n");
		if (delta_trigger >= 0.0)
		    delta = delta_next;
		return 1;
	    }
	}
#if !SCA
	if (DEBUG) {
	    // check if norm_g was computed correctly from g2
	    dgemv('T', m, n, Real(1.0), J.p, m, &r[0], 1, Real(0.0), &g[0], 1);
	    printf("        debugging:  norm_g:  %s  %s\n",
		   str(norm_g), str(norm_g - compute_norm(g)));
	}
#endif
	// update x, delta; compute new r
	int ret = check_reduction();
	if (DEBUG)
	    printf("ret = %d\n", ret);
	fflush(stdout);
	if (ret == 2) continue;
	return ret; // ret == 1 || ret == 0
    }
    if (verbose > 1)
	cout << "delta:  maxfev reached\n";
    fflush(stdout);
    return 0;
}


int levmar::check_reduction()    // updates x, delta
{
    // df = g*p + .5*(J*p)*(J*p) is known from compute_p()
    
    if (df >= 0) { // roundoff error can probably cause this
	if (verbose>1)
	    printf("df >= 0: delta %s, |p|: %s, df %s\n",
		   str(delta), str(norm_p2), str(df));
	return (stepsJ>0) ? 1 : 0; // 1: compute J,  0: done
    }

    old_x = x;
    update_x();

    if (x == old_x) {
	if (verbose > 1)
	    cout << "delta:  x == old_x\n";
	return (stepsJ>0) ? 1 : 0;
    }

    if (check_x() == 0) {
	x = old_x;
	if (verbose>1)
	    cout << "step rejected: check_x failed\n";
	fflush(stdout);

	if (stepsJ>0) {
	    if (delta_trigger >= 0.0)
		delta = delta_next;
	    return 1; // try again with new J
	}	
	delta = 0.375*norm_p2;
	return 2;  // J is correct. try again with smaller delta
    }

    old_f = f;
    old_r = r;
    old_norm_g = norm_g;

    compute_f(); nfev++;

#if SCA
    {
	int flags[4] = { 0 };
	flags[1] = 4;
	MPI_Bcast(flags, 4, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	scaM1->distribute(&r[0]);
	scaUU->apply_to_vec( *scaM1, *scaN1, 'T' );
	scaN1->consolidate(&g[0]);
	scale_vec(g, Sig);
	// no need to also apply VT since we only care about the norm here
    }
#else
    // note: the formula for g is r(x+p)*J(x0), where x0 is the
    //       value of x when stepsJ was zero.  (Jacobian is lagged)
    dgemv('T', m, n, Real(1.0), J.p, m, &r[0], 1, Real(0.0), &g[0], 1);
#endif

    norm_g = compute_norm(g);

    Real fmin = f;
    Real norm_g_min = norm_g;
    if (old_f < f) {
	fmin = old_f;
	norm_g_min = old_norm_g; // note: min on f, not g
    }
    Real rho = (f - old_f) / df; // df is negative
    
    // reset delta0, which is used in solve2 to pick new search directions
    if (update_delta0 && update_delta0_counter) {
	if (update_delta0_counter == 1)
	    delta0 = delta; // occurs 2nd time this code block is reached
	update_delta0_counter--;
    }

    // We could use norm_g0 instead of max(norm_g0,norm_g_min) if we
    // wanted to prevent roundoff regime from turning off due to
    // norm_g growing.  But it seems reasonable to toggle back to
    // normal if a larger gradient is encountered (even though the
    // Jacobian is lagged)

    if ( fmin < f_tol  ||  (abs(df) < df_tol*fmin && stepsJ==0)  ||
	 delta < g_tol ||  max(norm_g0,norm_g_min) < g_tol )
	return check_reduction_roundoff_regime(rho); // returns 0 or 2

    report_trust_reg(rho);

    if (rho > 0.0) { // accept the step
	if ( output_intermediate_successes )
	    dump(); // output successful iteration to file
	if ( ++stepsJ < max_stepsJ )
	    return update_delta(rho,2); // returns 1 or 2
	else
	    return update_delta(rho,1); // returns 1
    }

    x = old_x;
    f = old_f;
    r = old_r;
    norm_g = old_norm_g;
    // if (verbose>1) printf("step rejected\n");

    return update_delta(rho,2); // returns 1 or 2
}


int levmar::check_reduction_roundoff_regime(Real rho)
{
    report_trust_reg(rho, "roundoff regime");

    // do not recompute J (return 1) in roundoff regime
    if (f >= old_f) {
	x = old_x;
	f = old_f;
	r = old_r;
	norm_g = old_norm_g;
	if (verbose > 1)
	    printf("delta:  f >= old_f\n");
	return 0; // done
    }
    if (f > 0.999999999999*old_f) {
	// f < old_f, but not by much
	if (verbose > 1)
	    printf("delta:  f > .999*old_f\n");
	return 0; // done
    }
    if (f == 0.0) {
	if (verbose > 1)
	    printf("delta:  f == 0.0\n");
	return 0; // done
    }
    // f decreased significantly, so keep minimizing
    //  use_trigger=0 in roundoff regime, so update_delta will not return 1
    if ( ++stepsJ < max_stepsJ )
	return update_delta(rho,2,0); // returns 2
    else
	return update_delta(rho,0,0); // returns 0
}

int levmar::update_delta(Real rho, int ret0, int use_trigger)
{
    double rho0 = 0.250;  //  .9    .5   .25
    double rho1 = 0.850;  //  .975  .9   .75
    double   a0 = 0.375;  //  .5    .5   .25
    double   a1 = 1.875;  // 1.414 1.414  2
    double  np2 = 0.900;  //  .975  .9   .9
    double   a2 = 0.200;

    if (rho < rho0) {
	if (delta_trigger<0.0) {
	    // J was just computed if stepsJ==0 or stepsJ==1 && rho>0.0
	    // ~(A || B && C) = ~A && ( ~B || ~C) = ~A && ~B || ~A && ~C 
	    // if ( stepsJ>1  ||  (stepsJ>0 && rho<=0.0) ) {

	    // better idea: stepsJ>0 triggers immediately if a step is
	    //   accepted but the radius is rejected. otherwise works as before

	    if (stepsJ>0) {
		// delta = first rejected radius after an accepted step
		//   note: the radius was rejected, not necessarily the step
		//   note: max(a0^2, a0^3*a1^2) < a2 < min(a0, a0^2*a1)
		delta_trigger = delta*a2;
		// delta_next*a1 = first rejected radius after an accepted step
		delta_next = delta/a1;
	    }
	}
	delta = a0 * norm_p2;

	// use_trigger=0 in roundoff regime
	//   (allows delta to get really small without returning 1)
	if (use_trigger  &&  delta < delta_trigger) {
	    delta = delta_next;
	    return 1;
	}
    }
    if ( rho > rho1  &&  norm_p2/delta > np2 )
	delta = min(a1*delta, delta_max);

    return ret0;
}


void levmar::report_trust_reg(Real rho, const char *s) {
    if (verbose>1) {
	char step_buf[10] = "";
	char trigger_buf[30] = "";
	if (max_stepsJ>1) {
	    sprintf(step_buf, " %2d", stepsJ);
	    sprintf(trigger_buf, " trigger %8s,", str(delta_trigger,2));
	}
	printf("%4d%s, delta %9s,%s |p| %9s, lam %9s, df%10s, rho %9s, f%9s, "
	       "|g|%9s, newton%4d, smin%9s %s\n", num_steps, step_buf,
	       str(delta,3), trigger_buf,
	       str(norm_p2,3), str(lambda,3),
	       str(-df,4), str(rho,3), str(f,3), str(norm_g,3),
	       newton_iters, str(Sig[n-1],3), s);
    }
}


// fixme: probably better to use Sig*D instead of Lam. (avoid squaring)
void levmar::compute_p2()
{
    int i, i0;
    Real lam0(0);

    // check ordering of Lam
    for (i=0; i<n-1; i++)
	if (Lam[i]>Lam[i+1])
	    throw gen_err("error in order of Lam in levmar::compute_p2");

    // shift lambdas to avoid cancellation in lambda+Lam[0]
    if (Lam[0] < 0) {
	lam0 = Lam[0];
	Lam[0] = 0.0;
	for (i=1; i<n; i++)
	    Lam[i] -= lam0;
    }

    // find first index such that Lam[i0] != Lam[0]
    for (i0=1; i0<n; i0++)
	if (Lam[i0]>Lam[0]) break;

    Real g0sq = 0.0;
    for (int i=0; i<i0; i++)
	g0sq += g2[i]*g2[i];

    int done = 0;
    int k, kmax = 100;
    Real phi, tmp, tmp2;
    Real dlam, dphi, Lam0_lam;
    Real delta1 = 1.0/delta;
    //lambda = lam0;
    if (lambda < 0)
	lambda = 0.0;


    // phi is increasing and concave down for lambda > -Lam[0]
    //   (actually lambda > 0 after the shift)

    for (k=0; k<kmax; k++) {

	if (g0sq > 0 && lambda == -Lam[0]) {

	    phi = 0.0;
	    dphi = 1.0/sqrt(g0sq);

	} else {

	    Lam0_lam = Lam[0] + lambda;
	    if (Lam0_lam == 0)
		Lam0_lam = 1.0; // g0sq = 0 in this case
	    
	    phi = g0sq;
	    dphi = g0sq;
	    for (int i=i0; i<n; i++) {
		tmp = Lam0_lam/(Lam[i] + lambda);
		tmp2 = tmp*g2[i];
		tmp2 *= tmp2;
		phi += tmp2;
		dphi += tmp*tmp2;
	    }
	    if (phi == 0) break;
	    phi = 1.0/sqrt(phi);
	    dphi *= phi*phi*phi;
	    phi *= Lam0_lam;
	}

	dlam = (delta1 - phi)/dphi;
	if (verbose > 2)
	    printf("%s  %s  (lambda, dlam)\n", str(lambda), str(dlam));

	if ( lambda == 0  &&  dlam<=0 )
	    break;

	lambda += dlam;
	if (lambda < 0)
	    lambda = 0.0;

	if (done)
	    break;

	if ( abs(dlam) <= newton_tol*lambda ||
	     abs(delta1 - phi) <= newton_tol*delta1 )
	    done = 1; // do one more iteration

    } // end for k

    if (k>=kmax) {
	printf("**** error: Newton failed to converge.  "
	       "dlam = %s, lambda = %s ****\n",
	       str(dlam,4), str(lambda,4));
	printf("   i       Lam           g2\n");
	for (int i=0; i<n; i++)
	    printf("%4d %13s %13s\n", i, str(Lam[i]+lam0), str(g2[i]));
	throw gen_err("newton");
    }
    for (int i=0; i<n; i++) {
	if (i<i0 && g0sq == 0) p2[i] = 0.0;
	else p2[i] = -g2[i]/(lambda + Lam[i]);
    }
    norm_p2 = compute_norm(p2);

    if (norm_p2 > 1.0001*delta)
	throw gen_err("constraint error in compute_p2");

    if (lam0 < 0) {
	for (int i=0; i<n; i++)
	    Lam[i] += lam0;
    }
    df = 0.0;
    for (int i=0; i<n; i++)
	df += p2[i]*p2[i]*Lam[i];
    df = g2*p2 + .5*df;
    // printf("df = %s\n", str(df));

    newton_iters = k;
}
