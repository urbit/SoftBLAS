#include "softblas.h"

//  Construct a Givens rotation that zeroes b: returns c, s, and overwrites
//  a with r = +/-sqrt(a^2+b^2) and b with the reconstruction scalar z.
//  Quad-precision port of srotg using the f128M_* (pointer-output) API; no
//  call aliases an input pointer with its output pointer.
void qrotg(float128_t *a, float128_t *b, float128_t *c, float128_t *s, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float128_t ZERO = { SB_REAL128L_ZERO,   SB_REAL128U_ZERO };
    const float128_t ONE  = { SB_REAL128L_ONE,    SB_REAL128U_ONE };
    const float128_t NEG  = { SB_REAL128L_NEGONE, SB_REAL128U_NEGONE };

    float128_t A = *a, B = *b;
    float128_t absA = f128_abs(A), absB = f128_abs(B);
    float128_t roe   = f128M_gt(&absA, &absB) ? A : B;
    float128_t scale;
    f128M_add(&absA, &absB, &scale);
    float128_t r = ZERO, z = ZERO;

    if (f128M_eq(&scale, &ZERO)) {
        *c = ONE; *s = ZERO;
    } else {
        float128_t as, bs, as2, bs2, sumsq, sq;
        f128M_div(&A, &scale, &as);
        f128M_div(&B, &scale, &bs);
        f128M_mul(&as, &as, &as2);
        f128M_mul(&bs, &bs, &bs2);
        f128M_add(&as2, &bs2, &sumsq);
        f128M_sqrt(&sumsq, &sq);
        f128M_mul(&scale, &sq, &r);
        float128_t sgn = f128M_lt(&roe, &ZERO) ? NEG : ONE;
        float128_t rsigned;
        f128M_mul(&sgn, &r, &rsigned);   // r = sign(roe) * |r|
        r = rsigned;
        f128M_div(&A, &r, c);
        f128M_div(&B, &r, s);
        float128_t C = *c, S = *s;
        z = ONE;
        if (f128M_gt(&absA, &absB)) z = S;
        if (f128M_ge(&absB, &absA) && f128M_ne(&C, &ZERO)) {
            f128M_div(&ONE, &C, &z);
        }
    }
    *a = r;
    *b = z;
}
