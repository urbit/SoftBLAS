#include "softblas.h"

//  Construct a modified Givens transformation (reference SROTMG). Given the
//  diagonal scaling D1, D2, the first component X1 of the X-vector, and the
//  first component y1 of the Y-vector, it computes the flag P[0] and the H
//  matrix entries P[1..4] that zero the second component of the rotated
//  vector, rescaling so the scale factors stay in [1/gamsq, gamsq].
void srotmg(float32_t *D1, float32_t *D2, float32_t *X1, const float32_t y1, float32_t *P, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float32_t ZERO   = { SB_REAL32_ZERO };
    const float32_t ONE    = { SB_REAL32_ONE };
    const float32_t NEGONE = { SB_REAL32_NEGONE };
    const float32_t NEGTWO = { SB_REAL32_NEGTWO };
    const float32_t gam    = { 0x45800000 };  // 4096
    const float32_t gamsq  = { 0x4b800000 };  // 4096^2 = 2^24
    const float32_t rgam   = { 0x39800000 };  // 1/4096 = 2^-12
    const float32_t rgamsq = { 0x33800000 };  // 2^-24

    float32_t flag;
    float32_t d1 = *D1, d2 = *D2, x1 = *X1;
    float32_t h11 = ZERO, h21 = ZERO, h12 = ZERO, h22 = ZERO;

    if (f32_lt(d1, ZERO)) {
        flag = NEGONE;
        d1 = ZERO; d2 = ZERO; x1 = ZERO;
        h11 = ZERO; h21 = ZERO; h12 = ZERO; h22 = ZERO;
    } else {
        float32_t p2 = f32_mul(d2, y1);
        if (f32_eq(p2, ZERO)) {
            *P = NEGTWO;
            return;
        }
        float32_t p1 = f32_mul(d1, x1);
        float32_t q2 = f32_mul(p2, y1);
        float32_t q1 = f32_mul(p1, x1);

        if (f32_gt(f32_abs(q1), f32_abs(q2))) {
            h21 = f32_div(f32_mul(NEGONE, y1), x1);
            h12 = f32_div(p2, p1);
            float32_t u = f32_sub(ONE, f32_mul(h12, h21));
            if (f32_le(u, ZERO)) {       // singular -> zero result
                flag = NEGONE;
                d1 = ZERO; d2 = ZERO; x1 = ZERO;
                h11 = ZERO; h21 = ZERO; h12 = ZERO; h22 = ZERO;
                goto store;
            }
            flag = ZERO;
            d1 = f32_div(d1, u);
            d2 = f32_div(d2, u);
            x1 = f32_mul(x1, u);
        } else {
            if (f32_lt(q2, ZERO)) {
                flag = NEGONE;
                d1 = ZERO; d2 = ZERO; x1 = ZERO;
                h11 = ZERO; h21 = ZERO; h12 = ZERO; h22 = ZERO;
                goto store;
            }
            flag = ONE;
            h11 = f32_div(p1, p2);
            h22 = f32_div(x1, y1);
            float32_t u = f32_add(ONE, f32_mul(h11, h22));
            float32_t tmp = f32_div(d2, u);
            d2 = f32_div(d1, u);
            d1 = tmp;
            x1 = f32_mul(y1, u);
        }

        //  Rescale D1 into [1/gamsq, gamsq].
        while (f32_ne(d1, ZERO) && f32_le(d1, rgamsq)) {
            if (f32_eq(flag, ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
            else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
            d1  = f32_mul(d1, gamsq);
            x1  = f32_mul(x1, rgam);
            h11 = f32_mul(h11, rgam);
            h12 = f32_mul(h12, rgam);
        }
        while (f32_ge(d1, gamsq)) {
            if (f32_eq(flag, ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
            else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
            d1  = f32_mul(d1, rgamsq);
            x1  = f32_mul(x1, gam);
            h11 = f32_mul(h11, gam);
            h12 = f32_mul(h12, gam);
        }

        //  Rescale D2 into [1/gamsq, gamsq].
        while (f32_ne(d2, ZERO) && f32_le(f32_abs(d2), rgamsq)) {
            if (f32_eq(flag, ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
            else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
            d2  = f32_mul(d2, gamsq);
            h21 = f32_mul(h21, rgam);
            h22 = f32_mul(h22, rgam);
        }
        while (f32_ge(f32_abs(d2), gamsq)) {
            if (f32_eq(flag, ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
            else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
            d2  = f32_mul(d2, rgamsq);
            h21 = f32_mul(h21, gam);
            h22 = f32_mul(h22, gam);
        }
    }

store:
    *D1 = d1;
    *D2 = d2;
    *X1 = x1;
    P[0] = flag;
    if (f32_eq(flag, NEGONE)) { P[1] = h11; P[2] = h21; P[3] = h12; P[4] = h22; }
    else if (f32_eq(flag, ZERO)) { P[2] = h21; P[3] = h12; }
    else { P[1] = h11; P[4] = h22; }   // flag == 1
}
