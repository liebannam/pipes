#include "curves.h"

int main(int argc, char *argv[])
{
    int stage_sizes[] = {784, 400, 200, 100,  50,  25, 6,
			  25,  50, 100, 200, 400, 784};
    const char *stage_type[] = {
	"logistic", "logistic", "logistic",
	"logistic", "logistic", "linear",
	"logistic", "logistic", "logistic",
	"logistic", "logistic", "linear"};
    int nb = 1;
    int nQmax = 2;   // number of columns of history to remember
    int nk = 1;      // # of cols for apply_JT (0,1,nb)

    curves z(5000, 12, stage_sizes, stage_type, nb, nQmax, nk);
    
    z.load_images("bdata", "");

    // z.compute_f();
    // z.chkder(6.0e-6, 0);

    z.set_delta(1.0); // initial trust region size
    z.set_report_r_stride(0);
    z.solve();
    // z.output("results");

    return 0;
}
