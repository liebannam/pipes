#include <sstream>
#include <stdio.h>
#include "str_double.h"
using namespace std;

#define MX 50

const char *str(double x, int n) {
    static char buf[MX][32];
    static char fmt[8];
    static int k=MX-1;
    k++; if (k==MX) k=0;
    if (n==0 || n>17) n = 17;
    sprintf(fmt, "%%.%dg", n);
    sprintf(buf[k], fmt, x);
    return buf[k];
}

void str_to_real(const char *s, double &a) {
    sscanf(s, "%lf", &a);
}
