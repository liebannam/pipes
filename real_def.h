#ifndef __REAL_DEF_H
#define __REAL_DEF_H

#ifndef MPFR
#define MPFR 0
#endif

#ifndef DDDD
#define DDDD 0
#endif

#include "float.h"

template<class T> inline T real_eps(T a);

template<>
inline double real_eps(double a) {
    return DBL_EPSILON;
}

#if MPFR

#include "gmpfrxx.h"
typedef mpfr_class Real;
template<>
inline mpfr_class real_eps(mpfr_class a) {
    return mpfr_class(1.0) >> ( mpfr_class::get_dprec()-1 );
}

#elif DDDD

#include <qd/dd_real.h>
typedef dd_real Real;

template<>
inline dd_real real_eps(dd_real a) {
    return DBL_EPSILON*DBL_EPSILON;
}

#else

typedef double Real;

#endif

#endif
