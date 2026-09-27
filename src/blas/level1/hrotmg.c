#include "softblas.h"

//  Construct a modified Givens transformation (reference SROTMG). Given the
//  diagonal scaling D1, D2, the first component X1 of the X-vector, and the
//  first component y1 of the Y-vector, it computes the flag P[0] and the H
//  matrix entries P[1..4] that zero the second component of the rotated
//  vector, rescaling so the scale factors stay in [1/gamsq, gamsq].
void hrotmg(float16_t *D1, float16_t *D2, float16_t *X1, const float16_t y1, float16_t *P, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float16_t ZERO   = { SB_REAL16_ZERO };
    const float16_t ONE    = { SB_REAL16_ONE };
    const float16_t NEGONE = { SB_REAL16_NEGONE };
    const float16_t NEGTWO = { SB_REAL16_NEGTWO };
    const float16_t gam    = { 0x4c00 };  // 16   (reduced: float16 cannot hold 2^24)
    const float16_t gamsq  = { 0x5c00 };  // 256  = gam^2
    const float16_t rgam   = { 0x2c00 };  // 1/16
    const float16_t rgamsq = { 0x1c00 };  // 1/256
    const float16_t INF    = { 0x7c00 };  // +infinity

    float16_t flag;
    float16_t d1 = *D1, d2 = *D2, x1 = *X1;
    float16_t h11 = ZERO, h21 = ZERO, h12 = ZERO, h22 = ZERO;

    if (f16_lt(d1, ZERO)) {
        flag = NEGONE;
        d1 = ZERO; d2 = ZERO; x1 = ZERO;
        h11 = ZERO; h21 = ZERO; h12 = ZERO; h22 = ZERO;
    } else {
        float16_t p2 = f16_mul(d2, y1);
        if (f16_eq(p2, ZERO)) {
            *P = NEGTWO;
            return;
        }
        float16_t p1 = f16_mul(d1, x1);
        float16_t q2 = f16_mul(p2, y1);
        float16_t q1 = f16_mul(p1, x1);

        if (f16_gt(f16_abs(q1), f16_abs(q2))) {
            h21 = f16_div(f16_mul(NEGONE, y1), x1);
            h12 = f16_div(p2, p1);
            float16_t u = f16_sub(ONE, f16_mul(h12, h21));
            if (f16_le(u, ZERO)) {       // singular -> zero result
                flag = NEGONE;
                d1 = ZERO; d2 = ZERO; x1 = ZERO;
                h11 = ZERO; h21 = ZERO; h12 = ZERO; h22 = ZERO;
            } else {
                flag = ZERO;
                d1 = f16_div(d1, u);
                d2 = f16_div(d2, u);
                x1 = f16_mul(x1, u);
            }
        } else {
            if (f16_lt(q2, ZERO)) {      // singular -> zero result
                flag = NEGONE;
                d1 = ZERO; d2 = ZERO; x1 = ZERO;
                h11 = ZERO; h21 = ZERO; h12 = ZERO; h22 = ZERO;
            } else {
                flag = ONE;
                h11 = f16_div(p1, p2);
                h22 = f16_div(x1, y1);
                float16_t u = f16_add(ONE, f16_mul(h11, h22));
                float16_t tmp = f16_div(d2, u);
                d2 = f16_div(d1, u);
                d1 = tmp;
                x1 = f16_mul(y1, u);
            }
        }

        //  Rescale D1, D2 into [1/gamsq, gamsq]. A singular case above leaves
        //  d1 = d2 = 0, for which these loops are no-ops.
        //  Rescale D1 into [1/gamsq, gamsq]. The guard is 0 < d1, not d1 != 0:
        //  a negative d1 stays negative under the gamsq scaling and stays <=
        //  rgamsq all the way to -inf, so d1 != 0 would spin forever. (d1 < 0
        //  is rejected on entry, but the flag = 1 branch can still swap a
        //  negative d2 in when q2 = (d2*y1)*y1 underflows to -0 and so slips
        //  past the q2 < 0 test above.)
        while (f16_lt(ZERO, d1) && f16_le(d1, rgamsq)) {
            if (f16_eq(flag, ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
            else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
            d1  = f16_mul(d1, gamsq);
            x1  = f16_mul(x1, rgam);
            h11 = f16_mul(h11, rgam);
            h12 = f16_mul(h12, rgam);
        }
        while (f16_lt(f16_abs(d1), INF) && f16_ge(d1, gamsq)) {
            if (f16_eq(flag, ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
            else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
            d1  = f16_mul(d1, rgamsq);
            x1  = f16_mul(x1, gam);
            h11 = f16_mul(h11, gam);
            h12 = f16_mul(h12, gam);
        }

        //  Rescale D2 into [1/gamsq, gamsq].
        while (f16_ne(d2, ZERO) && f16_le(f16_abs(d2), rgamsq)) {
            if (f16_eq(flag, ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
            else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
            d2  = f16_mul(d2, gamsq);
            h21 = f16_mul(h21, rgam);
            h22 = f16_mul(h22, rgam);
        }
        while (f16_lt(f16_abs(d2), INF) && f16_ge(f16_abs(d2), gamsq)) {
            if (f16_eq(flag, ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
            else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
            d2  = f16_mul(d2, rgamsq);
            h21 = f16_mul(h21, gam);
            h22 = f16_mul(h22, gam);
        }
    }

    //  Canonicalize the values this call writes (flag is always one of
    //  -2/-1/0/1, never a NaN). The p2 == 0 early return above writes only
    //  P[0] and leaves D1/D2/X1 alone, so it needs nothing here.
    *D1 = nan_unify_h(d1);
    *D2 = nan_unify_h(d2);
    *X1 = nan_unify_h(x1);
    P[0] = flag;
    if (f16_eq(flag, NEGONE)) {
        P[1] = nan_unify_h(h11); P[2] = nan_unify_h(h21);
        P[3] = nan_unify_h(h12); P[4] = nan_unify_h(h22);
    }
    else if (f16_eq(flag, ZERO)) { P[2] = nan_unify_h(h21); P[3] = nan_unify_h(h12); }
    else { P[1] = nan_unify_h(h11); P[4] = nan_unify_h(h22); }   // flag == 1
}
