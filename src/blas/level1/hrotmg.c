#include "softblas.h"

void hrotmg(float16_t *D1, float16_t *D2, float16_t *X1, const float16_t y1, float16_t *P)
{
    float16_t flag = *P;
    float16_t d1=(*D1), d2=(*D2), x1=(*X1);
    float16_t h11, h21, h12, h22, tmp, u, p1, p2, q1, q2;
    static const float16_t gam=0x6c00;  // 4096.0
    static const float16_t gamsq = 0x6c00 * 0x6c00;
    static const float16_t rgam = SB_REAL16_ONE / 0x6c00;
    static const float16_t rgamsq = SB_REAL16_ONE /
                                        (0x6c00 * 0x6c00);

    if (f16_lt(d1, SB_REAL16_ZERO))
    {
        *P = SB_REAL16_NEGONE;
        *D1 = *D2 = *X1 = P[1] = P[2] = P[3] = P[4] = SB_REAL16_ZERO;
        return;
    }

    p2 = f16_mul(d2, y1);
    if (f16_eq(p2, SB_REAL16_ZERO))
    {
        *P = SB_REAL16_NEGTWO;
        return;
    }

    p1 = f16_mul(d1, x1);
    q2 = f16_mul(p2, y1);
    q1 = f16_mul(p1, x1);
    if (f16_gt(f16_abs(q1), f16_abs(q2)))
    {
        h21 = f16_div(f16_mul(SB_REAL16_NEGONE, y1), x1);
        h12 = f16_div(p2, p1);
        u = f16_sub(SB_REAL16_ONE, f16_mul(h12, h21));
        if (f16_gt(u, SB_REAL16_ZERO))
        {
            flag = SB_REAL16_ZERO;
            d1 = f16_div(d1, u);
            d2 = f16_div(d2, u);
            x1 = f16_mul(x1, u);
        }
    }
    else
    {
        if (f16_lt(q2, SB_REAL16_ZERO))
        {
            *P = SB_REAL16_NEGONE;
            *D1 = *D2 = *X1 = P[1] = P[2] = P[3] = P[4] = SB_REAL16_ZERO;
            return;
        }
        flag = SB_REAL16_ONE;
        h11 = f16_div(p1, p2);
        h22 = f16_div(x1, y1);
        u = f16_add(SB_REAL16_ONE, f16_mul(h11, h22));
        tmp = f16_div(d2, u);
        d2 = f16_div(d1, u);
        d1 = tmp;
        x1 = f16_mul(y1, u);

    }

    if (f16_le(d1, rgamsq))
    {
        if (f16_ne(d1, SB_REAL16_ZERO))
        {
            if (f16_eq(flag, SB_REAL16_ZERO)) { flag = SB_REAL16_NEGONE; h11 = h22 = SB_REAL16_ONE; }
            else { flag = h21 = SB_REAL16_NEGONE; h12 = SB_REAL16_ONE; }
            while(1)
            {
                d1 = f16_mul(d1, gamsq);
                x1 = f16_mul(x1, rgam);
                h11 = f16_mul(h11, rgam);
                h12 = f16_mul(h12, rgam);
                if (f16_gt(d1, rgamsq))
                    break;
                h21 = SB_REAL16_NEGONE;
                h12 = SB_REAL16_ONE;
            }
        }
    }
    else if (f16_ge(d1, gamsq))
    {
        if (f16_eq(flag, SB_REAL16_ZERO)) { flag = SB_REAL16_NEGONE; h11 = h22 = SB_REAL16_ONE; }
        else { flag = h21 = SB_REAL16_NEGONE; h12 = SB_REAL16_ONE; }
        while (1)
        {
          d1 = f16_mul(d1, rgamsq);
          x1 = f16_mul(x1, gam);
          h11 = f16_mul(h11, gam);
          h12 = f16_mul(h12, gam);
          if (f16_lt(d1, gamsq))
            break;
          h21 = SB_REAL16_NEGONE;
          h12 = SB_REAL16_ONE;
        }
    }
    
    tmp = f16_abs(d2);
    if (f16_le(tmp, rgamsq))
    {
        if (!f16_eq(d2, SB_REAL16_ZERO))
        {
            if (f16_eq(flag, SB_REAL16_ZERO)) { flag = SB_REAL16_NEGONE; h11 = h22 = SB_REAL16_ONE; }
            else { flag = h21 = SB_REAL16_NEGONE; h12 = SB_REAL16_ONE; }
            if (f16_ge(d2, SB_REAL16_ZERO))
            while (1)
            {
                d2 = f16_mul(d2, gamsq);
                h21 = f16_mul(h21, rgam);
                h22 = f16_mul(h22, rgam);
                if (f16_gt(d2, rgamsq))
                    break;
                h21 = SB_REAL16_NEGONE;
                h12 = SB_REAL16_ONE;
            }
            else /* d2 < SB_REAL16_ZERO */
            {
                tmp = f16_mul(SB_REAL16_NEGONE, rgamsq);
                while (1)
                {
                    d2 = f16_mul(d2, gamsq);
                    h21 = f16_mul(h21, rgam);
                    h22 = f16_mul(h22, rgam);
                    if (f16_lt(d2, tmp))
                        break;
                    h21 = SB_REAL16_NEGONE;
                    h12 = SB_REAL16_ONE;
                }
            }
        }
    }
    else if f16_ge(tmp, gamsq) {
        if (f16_eq(flag, SB_REAL16_ZERO)) { flag = SB_REAL16_NEGONE; h11 = h22 = SB_REAL16_ONE; }
        else { flag = h21 = SB_REAL16_NEGONE; h12 = SB_REAL16_ONE; }
        if (f16_ge(d2, SB_REAL16_ZERO))
        while (1)
        {
            d2 = f16_mul(d2, rgamsq);
            h21 = f16_mul(h21, gam);
            h22 = f16_mul(h22, gam);
            if (f16_lt(d2, gamsq))
                break;
            h21 = SB_REAL16_NEGONE;
            h12 = SB_REAL16_ONE;
        }
        else /* d2 < SB_REAL16_ZERO */
        {
            tmp = f16_mul(SB_REAL16_NEGONE, gamsq);
            while (1)
            {
                d2 = f16_mul(d2, rgamsq);
                h21 = f16_mul(h21, gam);
                h22 = f16_mul(h22, gam);
                if (f16_gt(d2, tmp))
                    break;
                h21 = SB_REAL16_NEGONE;
                h12 = SB_REAL16_ONE;
            }
        }
    }
   *D1 = d1;
   *D2 = d2;
   *X1 = x1;
   *P = flag;
   if (flag == SB_REAL16_NEGONE) { P[1] = h11; P[2] = h21; P[3] = h12; P[4] = h22; }
   else if (flag == SB_REAL16_ZERO) { P[2] = h21; P[3] = h12; }
   else if (flag == SB_REAL16_ONE) { P[1] = h11; P[4] = h22; }
}
