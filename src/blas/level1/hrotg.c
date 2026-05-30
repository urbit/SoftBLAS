#include "softblas.h"

//  Construct a Givens rotation that zeroes b: returns c, s, and overwrites
//  a with r = +/-sqrt(a^2+b^2) and b with the reconstruction scalar z.
//  (The SoftFloat comparison/abs macros do not parenthesize their argument,
//  so every operand here is first loaded into a plain local.)
void hrotg(float16_t *a, float16_t *b, float16_t *c, float16_t *s, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float16_t ZERO = { SB_REAL16_ZERO };
    const float16_t ONE  = { SB_REAL16_ONE };
    const float16_t NEG  = { SB_REAL16_NEGONE };

    float16_t A = *a, B = *b;
    float16_t absA = f16_abs(A), absB = f16_abs(B);
    float16_t roe   = f16_gt(absA, absB) ? A : B;
    float16_t scale = f16_add(absA, absB);
    float16_t r, z;

    if (f16_eq(scale, ZERO)) {
        *c = ONE; *s = ZERO; r = ZERO; z = ZERO;
    } else {
        float16_t as = f16_div(A, scale);
        float16_t bs = f16_div(B, scale);
        r = f16_mul(scale, f16_sqrt(f16_add(f16_mul(as, as), f16_mul(bs, bs))));
        r = f16_mul(f16_lt(roe, ZERO) ? NEG : ONE, r);   // sign(roe) * |r|
        *c = f16_div(A, r);
        *s = f16_div(B, r);
        float16_t C = *c, S = *s;
        z = ONE;
        if (f16_gt(absA, absB)) z = S;
        if (f16_ge(absB, absA) && f16_ne(C, ZERO)) z = f16_div(ONE, C);
    }
    *a = r;
    *b = z;
}
