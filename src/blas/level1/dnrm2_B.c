#include "softblas.h"

//  Euclidean norm via E. Anderson's "safe scaling" algorithm (ACM TOMS Alg.
//  978, 2017), as used by Reference-LAPACK/OpenBLAS. Three fixed-threshold
//  accumulators (small / medium / big) replace the per-element division of the
//  classic running-scale nrm2, while preserving overflow/underflow safety.
//  Constants are 2^k chosen from the float64 exponent range.
float64_t dnrm2_B(uint64_t N, const float64_t *X, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float64_t ZERO = { SB_REAL64_ZERO };
    const float64_t ONE  = { SB_REAL64_ONE };
    const float64_t tsml = { 0x2000000000000000 };  // 2^-511
    const float64_t tbig = { 0x5e50000000000000 };  // 2^486
    const float64_t ssml = { 0x6180000000000000 };  // 2^537
    const float64_t sbig = { 0x1e50000000000000 };  // 2^-538

    if (N < 1 || incX < 1) return nan_unify_d(ZERO);

    float64_t asml = ZERO, amed = ZERO, abig = ZERO;
    int notbig = 1;
    for (uint64_t k = 0, ix = 0; k < N; k++, ix += incX) {
        float64_t ax = f64_abs(X[ix]);
        if (f64_gt(ax, tbig)) {
            float64_t t = f64_mul(ax, sbig);
            abig = f64_add(abig, f64_mul(t, t));
            notbig = 0;
        } else if (f64_lt(ax, tsml)) {
            if (notbig) {
                float64_t t = f64_mul(ax, ssml);
                asml = f64_add(asml, f64_mul(t, t));
            }
        } else {
            amed = f64_add(amed, f64_mul(ax, ax));
        }
    }

    float64_t scl, sumsq;
    if (f64_gt(abig, ZERO)) {
        if (f64_gt(amed, ZERO) || f64_ne(amed, amed)) {
            abig = f64_add(abig, f64_mul(f64_mul(amed, sbig), sbig));
        }
        scl = f64_div(ONE, sbig);
        sumsq = abig;
    } else if (f64_gt(asml, ZERO)) {
        if (f64_gt(amed, ZERO) || f64_ne(amed, amed)) {
            float64_t am = f64_sqrt(amed);
            float64_t as = f64_div(f64_sqrt(asml), ssml);
            float64_t ymin = f64_gt(as, am) ? am : as;
            float64_t ymax = f64_gt(as, am) ? as : am;
            float64_t r = f64_div(ymin, ymax);
            scl = ONE;
            sumsq = f64_mul(f64_mul(ymax, ymax), f64_add(ONE, f64_mul(r, r)));
        } else {
            scl = f64_div(ONE, ssml);
            sumsq = asml;
        }
    } else {
        scl = ONE;
        sumsq = amed;
    }

    return nan_unify_d(f64_mul(scl, f64_sqrt(sumsq)));
}
