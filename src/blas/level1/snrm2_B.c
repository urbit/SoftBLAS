#include "softblas.h"

//  Euclidean norm via E. Anderson's "safe scaling" algorithm (ACM TOMS Alg.
//  978, 2017), as used by Reference-LAPACK/OpenBLAS. Three fixed-threshold
//  accumulators (small / medium / big) replace the per-element division of the
//  classic running-scale nrm2, while preserving overflow/underflow safety.
//  Constants are 2^k chosen from the float32 exponent range.
float32_t snrm2_B(uint64_t N, const float32_t *X, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float32_t ZERO = { SB_REAL32_ZERO };
    const float32_t ONE  = { SB_REAL32_ONE };
    const float32_t tsml = { 0x20000000 };  // 2^-63
    const float32_t tbig = { 0x59800000 };  // 2^52
    const float32_t ssml = { 0x65000000 };  // 2^75
    const float32_t sbig = { 0x19800000 };  // 2^-76

    if (N < 1 || incX < 1 || SB_STRIDE_OVERFLOWS(N, incX)) return nan_unify_s(ZERO);

    float32_t asml = ZERO, amed = ZERO, abig = ZERO;
    int notbig = 1;
    for (uint64_t k = 0, ix = 0; k < N; k++, ix += incX) {
        float32_t ax = f32_abs(X[ix]);
        if (f32_gt(ax, tbig)) {
            float32_t t = f32_mul(ax, sbig);
            abig = f32_add(abig, f32_mul(t, t));
            notbig = 0;
        } else if (f32_lt(ax, tsml)) {
            if (notbig) {
                float32_t t = f32_mul(ax, ssml);
                asml = f32_add(asml, f32_mul(t, t));
            }
        } else {
            amed = f32_add(amed, f32_mul(ax, ax));
        }
    }

    float32_t scl, sumsq;
    if (f32_gt(abig, ZERO)) {
        if (f32_gt(amed, ZERO) || f32_ne(amed, amed)) {
            abig = f32_add(abig, f32_mul(f32_mul(amed, sbig), sbig));
        }
        scl = f32_div(ONE, sbig);
        sumsq = abig;
    } else if (f32_gt(asml, ZERO)) {
        if (f32_gt(amed, ZERO) || f32_ne(amed, amed)) {
            float32_t am = f32_sqrt(amed);
            float32_t as = f32_div(f32_sqrt(asml), ssml);
            float32_t ymin = f32_gt(as, am) ? am : as;
            float32_t ymax = f32_gt(as, am) ? as : am;
            float32_t r = f32_div(ymin, ymax);
            scl = ONE;
            sumsq = f32_mul(f32_mul(ymax, ymax), f32_add(ONE, f32_mul(r, r)));
        } else {
            scl = f32_div(ONE, ssml);
            sumsq = asml;
        }
    } else {
        scl = ONE;
        sumsq = amed;
    }

    return nan_unify_s(f32_mul(scl, f32_sqrt(sumsq)));
}
