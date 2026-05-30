#include "softblas.h"

//  Construct a Givens rotation that zeroes b: returns c, s, and overwrites
//  a with r = +/-sqrt(a^2+b^2) and b with the reconstruction scalar z.
//  (The SoftFloat comparison/abs macros do not parenthesize their argument,
//  so every operand here is first loaded into a plain local.)
void srotg(float32_t *a, float32_t *b, float32_t *c, float32_t *s, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float32_t ZERO = { SB_REAL32_ZERO };
    const float32_t ONE  = { SB_REAL32_ONE };
    const float32_t NEG  = { SB_REAL32_NEGONE };

    float32_t A = *a, B = *b;
    float32_t absA = f32_abs(A), absB = f32_abs(B);
    float32_t roe   = f32_gt(absA, absB) ? A : B;
    float32_t scale = f32_add(absA, absB);
    float32_t r, z;

    if (f32_eq(scale, ZERO)) {
        *c = ONE; *s = ZERO; r = ZERO; z = ZERO;
    } else {
        float32_t as = f32_div(A, scale);
        float32_t bs = f32_div(B, scale);
        r = f32_mul(scale, f32_sqrt(f32_add(f32_mul(as, as), f32_mul(bs, bs))));
        r = f32_mul(f32_lt(roe, ZERO) ? NEG : ONE, r);   // sign(roe) * |r|
        *c = f32_div(A, r);
        *s = f32_div(B, r);
        float32_t C = *c, S = *s;
        z = ONE;
        if (f32_gt(absA, absB)) z = S;
        if (f32_ge(absB, absA) && f32_ne(C, ZERO)) z = f32_div(ONE, C);
    }
    *a = r;
    *b = z;
}
