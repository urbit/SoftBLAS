#include "softblas.h"

//  Construct a modified Givens transformation (reference SROTMG), quad
//  precision. Given the diagonal scaling D1, D2, the first component X1 of the
//  X-vector, and the first component y1 of the Y-vector, it computes the flag
//  P[0] and the H matrix entries P[1..4] that zero the second component of the
//  rotated vector, rescaling so the scale factors stay in [1/gamsq, gamsq].
//  Port of srotmg using the f128M_* (pointer-output) API; no call aliases an
//  input pointer with its output pointer.
void qrotmg(float128_t *D1, float128_t *D2, float128_t *X1, const float128_t y1, float128_t *P, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float128_t ZERO   = { SB_REAL128L_ZERO,   SB_REAL128U_ZERO };
    const float128_t ONE    = { SB_REAL128L_ONE,    SB_REAL128U_ONE };
    const float128_t NEGONE = { SB_REAL128L_NEGONE, SB_REAL128U_NEGONE };
    const float128_t NEGTWO = { SB_REAL128L_NEGTWO, SB_REAL128U_NEGTWO };
    const float128_t gam    = { 0x0000000000000000, 0x400b000000000000 };  // 4096 = 2^12
    const float128_t gamsq  = { 0x0000000000000000, 0x4017000000000000 };  // 4096^2 = 2^24
    const float128_t rgam   = { 0x0000000000000000, 0x3ff3000000000000 };  // 1/4096 = 2^-12
    const float128_t rgamsq = { 0x0000000000000000, 0x3fe7000000000000 };  // 2^-24
    const float128_t INF    = { 0x0000000000000000, 0x7fff000000000000 };  // +infinity

    float128_t flag;
    float128_t d1 = *D1, d2 = *D2, x1 = *X1;
    float128_t h11 = ZERO, h21 = ZERO, h12 = ZERO, h22 = ZERO;

    if (f128M_lt(&d1, &ZERO)) {
        flag = NEGONE;
        d1 = ZERO; d2 = ZERO; x1 = ZERO;
    } else {
        float128_t p2;
        f128M_mul(&d2, &y1, &p2);
        if (f128M_eq(&p2, &ZERO)) {
            *P = NEGTWO;
            return;
        }
        float128_t p1, q2, q1;
        f128M_mul(&d1, &x1, &p1);
        f128M_mul(&p2, &y1, &q2);
        f128M_mul(&p1, &x1, &q1);

        float128_t aq1 = f128_abs(q1), aq2 = f128_abs(q2);
        if (f128M_gt(&aq1, &aq2)) {
            float128_t ny1, t, u;
            f128M_mul(&NEGONE, &y1, &ny1);
            f128M_div(&ny1, &x1, &h21);   // h21 = -y1/x1
            f128M_div(&p2, &p1, &h12);    // h12 = p2/p1
            f128M_mul(&h12, &h21, &t);
            f128M_sub(&ONE, &t, &u);      // u = 1 - h12*h21
            if (f128M_le(&u, &ZERO)) {    // singular -> zero result
                flag = NEGONE;
                d1 = ZERO; d2 = ZERO; x1 = ZERO;
                h11 = ZERO; h21 = ZERO; h12 = ZERO; h22 = ZERO;
            } else {
                flag = ZERO;
                float128_t t1, t2, t3;
                f128M_div(&d1, &u, &t1); d1 = t1;
                f128M_div(&d2, &u, &t2); d2 = t2;
                f128M_mul(&x1, &u, &t3); x1 = t3;
            }
        } else {
            if (f128M_lt(&q2, &ZERO)) {   // singular -> zero result
                flag = NEGONE;
                d1 = ZERO; d2 = ZERO; x1 = ZERO;
                h11 = ZERO; h21 = ZERO; h12 = ZERO; h22 = ZERO;
            } else {
                flag = ONE;
                f128M_div(&p1, &p2, &h11);  // h11 = p1/p2
                f128M_div(&x1, &y1, &h22);  // h22 = x1/y1
                float128_t t, u, tmp, nd2;
                f128M_mul(&h11, &h22, &t);
                f128M_add(&ONE, &t, &u);    // u = 1 + h11*h22
                f128M_div(&d2, &u, &tmp);   // tmp = d2/u
                f128M_div(&d1, &u, &nd2);   // d2  = d1/u
                d2 = nd2;
                d1 = tmp;
                f128M_mul(&y1, &u, &x1);    // x1 = y1*u
            }
        }

        //  Rescale D1, D2 into [1/gamsq, gamsq]. A singular case above leaves
        //  d1 = d2 = 0, for which these loops are no-ops. The _ge loops are
        //  guarded with a finiteness check: f128M_ge(d, gamsq) expands to
        //  !(f128M_lt(d, gamsq)), which is true for any non-finite d, so an
        //  Inf/NaN scale factor would otherwise spin them forever.
        {
            float128_t t;
            while (f128M_ne(&d1, &ZERO) && f128M_le(&d1, &rgamsq)) {
                if (f128M_eq(&flag, &ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
                else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
                f128M_mul(&d1, &gamsq, &t); d1 = t;
                f128M_mul(&x1, &rgam, &t);  x1 = t;
                f128M_mul(&h11, &rgam, &t); h11 = t;
                f128M_mul(&h12, &rgam, &t); h12 = t;
            }
            float128_t ad1 = f128_abs(d1);
            while (f128M_lt(&ad1, &INF) && f128M_ge(&d1, &gamsq)) {
                if (f128M_eq(&flag, &ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
                else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
                f128M_mul(&d1, &rgamsq, &t); d1 = t;
                f128M_mul(&x1, &gam, &t);    x1 = t;
                f128M_mul(&h11, &gam, &t);   h11 = t;
                f128M_mul(&h12, &gam, &t);   h12 = t;
                ad1 = f128_abs(d1);
            }

            float128_t ad2 = f128_abs(d2);
            while (f128M_ne(&d2, &ZERO) && f128M_le(&ad2, &rgamsq)) {
                if (f128M_eq(&flag, &ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
                else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
                f128M_mul(&d2, &gamsq, &t); d2 = t;
                f128M_mul(&h21, &rgam, &t); h21 = t;
                f128M_mul(&h22, &rgam, &t); h22 = t;
                ad2 = f128_abs(d2);
            }
            while (f128M_lt(&ad2, &INF) && f128M_ge(&ad2, &gamsq)) {
                if (f128M_eq(&flag, &ZERO)) { flag = NEGONE; h11 = ONE; h22 = ONE; }
                else { flag = NEGONE; h21 = NEGONE; h12 = ONE; }
                f128M_mul(&d2, &rgamsq, &t); d2 = t;
                f128M_mul(&h21, &gam, &t);   h21 = t;
                f128M_mul(&h22, &gam, &t);   h22 = t;
                ad2 = f128_abs(d2);
            }
        }
    }

    *D1 = d1;
    *D2 = d2;
    *X1 = x1;
    P[0] = flag;
    if (f128M_eq(&flag, &NEGONE)) { P[1] = h11; P[2] = h21; P[3] = h12; P[4] = h22; }
    else if (f128M_eq(&flag, &ZERO)) { P[2] = h21; P[3] = h12; }
    else { P[1] = h11; P[4] = h22; }   // flag == 1
}
