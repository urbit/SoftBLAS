#include "softblas.h"

//  Euclidean norm via E. Anderson's "safe scaling" algorithm (ACM TOMS Alg.
//  978, 2017). Quad-precision variant using the f128M_ pointer API. Three
//  fixed-threshold accumulators (small / medium / big) replace the per-element
//  division of the classic running-scale nrm2. Constants are 2^k from the
//  float128 exponent range.
float128_t qnrm2_B(uint64_t N, const float128_t *X, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float128_t ZERO = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };
    const float128_t ONE  = { SB_REAL128L_ONE,  SB_REAL128U_ONE };
    const float128_t tsml = { 0x0000000000000000, 0x2000000000000000 };  // 2^-8191
    const float128_t tbig = { 0x0000000000000000, 0x5fc7000000000000 };  // 2^8136
    const float128_t ssml = { 0x0000000000000000, 0x6036000000000000 };  // 2^8247
    const float128_t sbig = { 0x0000000000000000, 0x1fc7000000000000 };  // 2^-8248

    float128_t norm = ZERO;
    if (N < 1 || incX < 1) { nan_unify_q(&norm); return norm; }

    float128_t asml = ZERO, amed = ZERO, abig = ZERO;
    float128_t ax, t, tmp;
    int notbig = 1;
    for (uint64_t k = 0, ix = 0; k < N; k++, ix += incX) {
        ax = f128_abs(X[ix]);
        if (f128M_gt(&ax, &tbig)) {
            f128M_mul(&ax, &sbig, &t);
            f128M_mul(&t, &t, &tmp);
            f128M_add(&abig, &tmp, &abig);
            notbig = 0;
        } else if (f128M_lt(&ax, &tsml)) {
            if (notbig) {
                f128M_mul(&ax, &ssml, &t);
                f128M_mul(&t, &t, &tmp);
                f128M_add(&asml, &tmp, &asml);
            }
        } else {
            f128M_mul(&ax, &ax, &tmp);
            f128M_add(&amed, &tmp, &amed);
        }
    }

    float128_t scl, sumsq;
    if (f128M_gt(&abig, &ZERO)) {
        if (f128M_gt(&amed, &ZERO) || f128M_ne(&amed, &amed)) {
            f128M_mul(&amed, &sbig, &t);
            f128M_mul(&t, &sbig, &tmp);
            f128M_add(&abig, &tmp, &abig);
        }
        f128M_div(&ONE, &sbig, &scl);
        sumsq = abig;
    } else if (f128M_gt(&asml, &ZERO)) {
        if (f128M_gt(&amed, &ZERO) || f128M_ne(&amed, &amed)) {
            float128_t am, as, ymin, ymax, r;
            f128M_sqrt(&amed, &am);
            f128M_sqrt(&asml, &t);
            f128M_div(&t, &ssml, &as);
            if (f128M_gt(&as, &am)) { ymin = am; ymax = as; }
            else                   { ymin = as; ymax = am; }
            f128M_div(&ymin, &ymax, &r);
            f128M_mul(&r, &r, &tmp);          // r^2
            f128M_add(&ONE, &tmp, &t);        // 1 + r^2
            f128M_mul(&ymax, &ymax, &tmp);    // ymax^2
            f128M_mul(&tmp, &t, &sumsq);      // ymax^2 * (1 + r^2)
            scl = ONE;
        } else {
            f128M_div(&ONE, &ssml, &scl);
            sumsq = asml;
        }
    } else {
        scl = ONE;
        sumsq = amed;
    }

    f128M_sqrt(&sumsq, &t);
    f128M_mul(&scl, &t, &norm);
    nan_unify_q(&norm);
    return norm;
}
