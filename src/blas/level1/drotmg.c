#include "softblas.h"

//  Construct a modified Givens transformation (reference SROTMG). Given the
//  diagonal scaling D1, D2, the first component X1 of the X-vector, and the
//  first component y1 of the Y-vector, it computes the flag P[0] and the H
//  matrix entries P[1..4] that zero the second component of the rotated
//  vector, rescaling so the scale factors stay in [1/gamsq, gamsq].
void drotmg(float64_t *D1, float64_t *D2, float64_t *X1, const float64_t y1, float64_t *P, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float64_t ZERO   = { SB_REAL64_ZERO };
    const float64_t ONE    = { SB_REAL64_ONE };
    const float64_t NEGONE = { SB_REAL64_NEGONE };
    const float64_t NEGTWO = { SB_REAL64_NEGTWO };
    const float64_t gam    = { 0x40b0000000000000 };  // 4096
    const float64_t gamsq  = { 0x4170000000000000 };  // 4096^2 = 2^24
    const float64_t rgam   = { 0x3f30000000000000 };  // 1/4096 = 2^-12
    const float64_t rgamsq = { 0x3e70000000000000 };  // 2^-24

    float64_t flag;
    float64_t d1 = *D1, d2 = *D2, x1 = *X1;
    float64_t h11 = ZERO, h21 = ZERO, h12 = ZERO, h22 = ZERO;

    if (f64_lt(d1, ZERO)) {
        flag = NEGONE;
        d1 = ZERO; d2 = ZERO; x1 = ZERO;
        h11 = ZERO; h21 = ZERO; h12 = ZERO; h22 = ZERO;
    } else {
        float64_t p2 = f64_mul(d2, y1);
        if (f64_eq(p2, ZERO)) {
            *P = NEGTWO;
            return;
        }
        float64_t p1 = f64_mul(d1, x1);
        float64_t q2 = f64_mul(p2, y1);
        float64_t q1 = f64_mul(p1, x1);

        if (f64_gt(f64_abs(q1), f64_abs(q2))) {
            h21 = f64_div(f64_mul(NEGONE, y1), x1);
            h12 = f64_div(p2, p1);
            float64_t u = f64_sub(ONE, f64_mul(h12, h21));
            if (f64_le(u, ZERO)) {       // singular -> zero result
                flag = NEGONE;
                d1 = ZERO; d2 = ZERO; x1 = ZERO;
                h11 = ZERO; h21 = ZERO; h12 = ZERO; h22 = ZERO;
                goto store;
            }
            flag = ZERO;
            d1 = f64_div(d1, u);
            d2 = f64_div(d2, u);
            x1 = f64_mul(x1, u);
        } else {
            if (f64_lt(q2, ZERO)) {
                flag = NEGONE;
                d1 = ZERO; d2 = ZERO; x1 = ZERO;
                h11 = ZERO; h21 = ZERO; h12 = ZERO; h22 = ZERO;
                goto store;
            }
            flag = ONE;
            h11 = f64_div(p1, p2);
            h22 = f64_div(x1, y1);
            float64_t u = f64_add(ONE, f64_mul(h11, h22));
            float64_t tmp = f64_div(d2, u);
            d2 = f64_div(d1, u);
            d1 = tmp;
            x1 = f64_mul(y1, u);
        }

        //  Rescale D1 into [1/gamsq, gamsq].
        while (f64_ne(d1, ZERO) && f64_le(d1, rgamsq)) {
            if (f64_eq(flag, ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
            else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
            d1  = f64_mul(d1, gamsq);
            x1  = f64_mul(x1, rgam);
            h11 = f64_mul(h11, rgam);
            h12 = f64_mul(h12, rgam);
        }
        while (f64_ge(d1, gamsq)) {
            if (f64_eq(flag, ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
            else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
            d1  = f64_mul(d1, rgamsq);
            x1  = f64_mul(x1, gam);
            h11 = f64_mul(h11, gam);
            h12 = f64_mul(h12, gam);
        }

        //  Rescale D2 into [1/gamsq, gamsq].
        while (f64_ne(d2, ZERO) && f64_le(f64_abs(d2), rgamsq)) {
            if (f64_eq(flag, ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
            else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
            d2  = f64_mul(d2, gamsq);
            h21 = f64_mul(h21, rgam);
            h22 = f64_mul(h22, rgam);
        }
        while (f64_ge(f64_abs(d2), gamsq)) {
            if (f64_eq(flag, ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
            else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
            d2  = f64_mul(d2, rgamsq);
            h21 = f64_mul(h21, gam);
            h22 = f64_mul(h22, gam);
        }
    }

store:
    *D1 = d1;
    *D2 = d2;
    *X1 = x1;
    P[0] = flag;
    if (f64_eq(flag, NEGONE)) { P[1] = h11; P[2] = h21; P[3] = h12; P[4] = h22; }
    else if (f64_eq(flag, ZERO)) { P[2] = h21; P[3] = h12; }
    else { P[1] = h11; P[4] = h22; }   // flag == 1
}
