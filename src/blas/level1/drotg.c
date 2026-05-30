#include "softblas.h"

//  Construct a Givens rotation that zeroes b: returns c, s, and overwrites
//  a with r = +/-sqrt(a^2+b^2) and b with the reconstruction scalar z.
//  (The SoftFloat comparison/abs macros do not parenthesize their argument,
//  so every operand here is first loaded into a plain local.)
void drotg(float64_t *a, float64_t *b, float64_t *c, float64_t *s, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float64_t ZERO = { SB_REAL64_ZERO };
    const float64_t ONE  = { SB_REAL64_ONE };
    const float64_t NEG  = { SB_REAL64_NEGONE };

    float64_t A = *a, B = *b;
    float64_t absA = f64_abs(A), absB = f64_abs(B);
    float64_t roe   = f64_gt(absA, absB) ? A : B;
    float64_t scale = f64_add(absA, absB);
    float64_t r, z;

    if (f64_eq(scale, ZERO)) {
        *c = ONE; *s = ZERO; r = ZERO; z = ZERO;
    } else {
        float64_t as = f64_div(A, scale);
        float64_t bs = f64_div(B, scale);
        r = f64_mul(scale, f64_sqrt(f64_add(f64_mul(as, as), f64_mul(bs, bs))));
        r = f64_mul(f64_lt(roe, ZERO) ? NEG : ONE, r);   // sign(roe) * |r|
        *c = f64_div(A, r);
        *s = f64_div(B, r);
        float64_t C = *c, S = *s;
        z = ONE;
        if (f64_gt(absA, absB)) z = S;
        if (f64_ge(absB, absA) && f64_ne(C, ZERO)) z = f64_div(ONE, C);
    }
    *a = r;
    *b = z;
}
