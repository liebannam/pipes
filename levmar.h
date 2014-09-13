#ifndef __LEVMAR
#define __LEVMAR

#include "real_def.h"
#include "mp_mat.h"
#include <vector>
using namespace std;

#if DDDD
#include "str_dd.h"
#include "lapack_dd.h"
#else
#include "str_double.h"
#include "lapack.h"
#endif

// be careful. make sure all files that include sca.h
//   are consistent in their definition of SCA
#if SCA
#include "sca.h"
#endif

class levmar {
public:
    int m;
    int n;

    Real f;
    Real df;

    vector<Real> x;  // length n
    vector<Real> p;  // length n
    vector<Real> r;  // length m
    vector<Real> g;  // length n

    vector<Real> old_x;
    vector<Real> old_r;
    Real old_f;
    Real old_norm_g;

    Real gammaD;     // D[i] = Sig[i]^(-gammaD), 0.0 <= gammaD < 1.0
    vector<Real> D;  // length n   trust region scaling weights
    vector<Real> p2; // length n   p2 = svd coordinates of p, p2 = D^{-1}V'*p
    vector<Real> g2; // length n   g2 = svd coordinates of g, g2 = DV'*g
    vector<Real> tmp;

    Real norm_g0; // = norm_g on first pass of take_step loop (when stepsJ==0)
    Real norm_g; // g' = r'*J
    Real norm_p2;
    Real norm_r;  // user should *not* update norm_r in compute_f

    mp_mat<Real> J;    // m by n

    mp_mat<Real>  UU;  // m by n
    mp_mat<Real>  VT;  // n by n
    vector<Real>  Sig; // n
    vector<Real>  Lam; // n

    Real lambda;
    int verbose;

    bool update_delta0; // sets delta0 to delta on 2nd call of check_reduction()
    int update_delta0_counter; // helper for update_delta0
    Real delta;
    Real delta0;
    Real delta_max;
    Real delta_trigger;
    Real delta_next;

    int nfev;
    int nJev;
    int maxfev;

    int num_steps;
    int stepsJ;
    int max_stepsJ;
    int report_r_stride;
    int newton_iters;

    Real f_tol;
    Real g_tol;
    Real df_tol;
    Real newton_tol;
    Real diag_tol;
    Real sig_tol;

#if SCA
    dist_mat *scaJq;
    dist_mat *scaUU;
    dist_mat *scaVT;
    dist_mat *scaN1; // n by 1
    dist_mat *scaN2; // n by 1
    dist_mat *scaM1; // m by 1
    void set_dist_mat_ptrs(dist_mat *a, dist_mat *b, dist_mat *c,
			   dist_mat *d, dist_mat *e, dist_mat *f);
#endif

    levmar(int m, int n);

    virtual ~levmar() { }

    void set_x(Real *xx);
    void set_gammaD(Real gD);
    void set_init_params();

    void update_x(); // add p to x
    virtual int check_x() { return 1; }
    void compute_diag(int nD); // D[i] = Sig[i]^(-gammaD), 0.0 <= gammaD < 1.0

    void chkder(Real dx = 6.0e-6);
    void solve(int skipJ=0); // 1 means J already set up on first pass

    void compute_f(); // 0.5*r'*r
    void compute_g(); // J'*r

    virtual void compute_r() = 0;
    virtual void compute_J() = 0;

    // return 0 if x==old_x or a step is rejected in roundoff regime
    //               or stepsJ reaches max_stepsJ in roundoff regime
    // return 1 if not in roundoff regime and:
    //          max_stepsJ is reached, or if a step has been
    //          accepted and, later, delta drops below delta_trigger,
    //          or if stepsJ>0, and df is tiny
    //          (1 means the next batch of J should be computed)
    // return 2 if a step is rejected on a freshly computed Jacobian,
    //          or if a step is accepted or rejected without triggering 1
    int take_step();

    int check_reduction(); // update x, delta;  0:done, 1:update J, 2:re-use J
    int check_reduction_roundoff_regime(Real rho);
    int update_delta(Real rho, int ret0, int use_trigger=1);
    void report_trust_reg(Real rho, const char *s = "");
    void compute_p2();

    // void load(const char *fname); // fixme
    virtual void dump();
    void report_r(FILE *fp = stdout);

    // 0 by default, calls dump() each time a step is accepted
    int output_intermediate_successes;

    void set_maxfev(int mfev) { maxfev = mfev; }
    void set_delta_max(Real dmax) { delta_max = dmax; }
    // should call set_delta or reset_delta between calls of solve()
    void set_delta(Real d) { delta = d; delta0 = d; }
    void reset_delta() { delta = delta0; }
    void set_verbose(int vb) { verbose = vb; }
    void set_newton_tol(Real nt) { newton_tol = nt; }
    void set_f_tol(Real ft) { f_tol = ft; }
    void set_g_tol(Real gt) { g_tol = gt; }
    void set_df_tol(Real dft) { df_tol = dft; }
    void set_max_stepsJ(int k) { max_stepsJ = k; }
};

#endif
