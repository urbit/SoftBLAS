#include "softblas.h"

//  Euclidean norm via E. Anderson's "safe scaling" algorithm (ACM TOMS Alg.
//  978, 2017), as used by Reference-LAPACK/OpenBLAS. Three fixed-threshold
//  accumulators (small / medium / big) replace the per-element division of the
//  classic running-scale nrm2, while preserving overflow/underflow safety.
//  Constants are 2^k chosen from the float16 exponent range.
float16_t hnrm2_B(uint64_t N, const float16_t *X, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float16_t ZERO = { SB_REAL16_ZERO };
    const float16_t ONE  = { SB_REAL16_ONE };
    const float16_t tsml = { 0x2000 };  // 2^-7
    const float16_t tbig = { 0x4800 };  // 2^3
    const float16_t ssml = { 0x6c00 };  // 2^12
    const float16_t sbig = { 0x0800 };  // 2^-13

    if (N < 1 || incX < 1) return nan_unify_h(ZERO);

    float16_t asml = ZERO, amed = ZERO, abig = ZERO;
    int notbig = 1;
    for (uint64_t k = 0, ix = 0; k < N; k++, ix += incX) {
        float16_t ax = f16_abs(X[ix]);
        if (f16_gt(ax, tbig)) {
            float16_t t = f16_mul(ax, sbig);
            abig = f16_add(abig, f16_mul(t, t));
            notbig = 0;
        } else if (f16_lt(ax, tsml)) {
            if (notbig) {
                float16_t t = f16_mul(ax, ssml);
                asml = f16_add(asml, f16_mul(t, t));
            }
        } else {
            amed = f16_add(amed, f16_mul(ax, ax));
        }
    }

    float16_t scl, sumsq;
    if (f16_gt(abig, ZERO)) {
        if (f16_gt(amed, ZERO) || f16_ne(amed, amed)) {
            abig = f16_add(abig, f16_mul(f16_mul(amed, sbig), sbig));
        }
        scl = f16_div(ONE, sbig);
        sumsq = abig;
    } else if (f16_gt(asml, ZERO)) {
        if (f16_gt(amed, ZERO) || f16_ne(amed, amed)) {
            float16_t am = f16_sqrt(amed);
            float16_t as = f16_div(f16_sqrt(asml), ssml);
            float16_t ymin = f16_gt(as, am) ? am : as;
            float16_t ymax = f16_gt(as, am) ? as : am;
            float16_t r = f16_div(ymin, ymax);
            scl = ONE;
            sumsq = f16_mul(f16_mul(ymax, ymax), f16_add(ONE, f16_mul(r, r)));
        } else {
            scl = f16_div(ONE, ssml);
            sumsq = asml;
        }
    } else {
        scl = ONE;
        sumsq = amed;
    }

    return nan_unify_h(f16_mul(scl, f16_sqrt(sumsq)));
}
