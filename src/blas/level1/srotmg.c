#include "softblas.h"

void srotmg(float32_t *D1, float32_t *D2, float32_t *X1, const float32_t y1, float32_t *P)
{
    float32_t flag = *P;
    float32_t d1=(*D1), d2=(*D2), x1=(*X1);
    float32_t h11, h21, h12, h22, tmp, u, p1, p2, q1, q2;
    static const float32_t gam=0x45800000;  // 4096.0
    static const float32_t gamsq = 0x45800000 * 0x45800000;
    static const float32_t rgam = SB_REAL32_ONE / 0x45800000;
    static const float32_t rgamsq = SB_REAL32_ONE /
                                        (0x45800000 * 0x45800000);

    if (f32_lt(d1, SB_REAL32_ZERO))
    {
        *P = SB_REAL32_NEGONE;
        *D1 = *D2 = *X1 = P[1] = P[2] = P[3] = P[4] = SB_REAL32_ZERO;
        return;
    }

    p2 = f32_mul(d2, y1);
    if (f32_eq(p2, SB_REAL32_ZERO))
    {
        *P = SB_REAL32_NEGTWO;
        return;
    }

    p1 = f32_mul(d1, x1);
    q2 = f32_mul(p2, y1);
    q1 = f32_mul(p1, x1);
    if (f32_gt(f32_abs(q1), f32_abs(q2)))
    {
        h21 = f32_div(f32_mul(SB_REAL32_NEGONE, y1), x1);
        h12 = f32_div(p2, p1);
        u = f32_sub(SB_REAL32_ONE, f32_mul(h12, h21));
        if (f32_gt(u, SB_REAL32_ZERO))
        {
            flag = SB_REAL32_ZERO;
            d1 = f32_div(d1, u);
            d2 = f32_div(d2, u);
            x1 = f32_mul(x1, u);
        }
    }
    else
    {
        if (f32_lt(q2, SB_REAL32_ZERO))
        {
            *P = SB_REAL32_NEGONE;
            *D1 = *D2 = *X1 = P[1] = P[2] = P[3] = P[4] = SB_REAL32_ZERO;
            return;
        }
        flag = SB_REAL32_ONE;
        h11 = f32_div(p1, p2);
        h22 = f32_div(x1, y1);
        u = f32_add(SB_REAL32_ONE, f32_mul(h11, h22));
        tmp = f32_div(d2, u);
        d2 = f32_div(d1, u);
        d1 = tmp;
        x1 = f32_mul(y1, u);

    }

    if (f32_le(d1, rgamsq))
    {
        if (f32_ne(d1, SB_REAL32_ZERO))
        {
            if (f32_eq(flag, SB_REAL32_ZERO)) { flag = SB_REAL32_NEGONE; h11 = h22 = SB_REAL32_ONE; }
            else { flag = h21 = SB_REAL32_NEGONE; h12 = SB_REAL32_ONE; }
            while(1)
            {
                d1 = f32_mul(d1, gamsq);
                x1 = f32_mul(x1, rgam);
                h11 = f32_mul(h11, rgam);
                h12 = f32_mul(h12, rgam);
                if (f32_gt(d1, rgamsq))
                    break;
                h21 = SB_REAL32_NEGONE;
                h12 = SB_REAL32_ONE;
            }
        }
    }
    else if (f32_ge(d1, gamsq))
    {
        if (f32_eq(flag, SB_REAL32_ZERO)) { flag = SB_REAL32_NEGONE; h11 = h22 = SB_REAL32_ONE; }
        else { flag = h21 = SB_REAL32_NEGONE; h12 = SB_REAL32_ONE; }
        while (1)
        {
          d1 = f32_mul(d1, rgamsq);
          x1 = f32_mul(x1, gam);
          h11 = f32_mul(h11, gam);
          h12 = f32_mul(h12, gam);
          if (f32_lt(d1, gamsq))
            break;
          h21 = SB_REAL32_NEGONE;
          h12 = SB_REAL32_ONE;
        }
    }
    
    tmp = f32_abs(d2);
    if (f32_le(tmp, rgamsq))
    {
        if (!f32_eq(d2, SB_REAL32_ZERO))
        {
            if (f32_eq(flag, SB_REAL32_ZERO)) { flag = SB_REAL32_NEGONE; h11 = h22 = SB_REAL32_ONE; }
            else { flag = h21 = SB_REAL32_NEGONE; h12 = SB_REAL32_ONE; }
            if (f32_ge(d2, SB_REAL32_ZERO))
            while (1)
            {
                d2 = f32_mul(d2, gamsq);
                h21 = f32_mul(h21, rgam);
                h22 = f32_mul(h22, rgam);
                if (f32_gt(d2, rgamsq))
                    break;
                h21 = SB_REAL32_NEGONE;
                h12 = SB_REAL32_ONE;
            }
            else /* d2 < SB_REAL32_ZERO */
            {
                tmp = f32_mul(SB_REAL32_NEGONE, rgamsq);
                while (1)
                {
                    d2 = f32_mul(d2, gamsq);
                    h21 = f32_mul(h21, rgam);
                    h22 = f32_mul(h22, rgam);
                    if (f32_lt(d2, tmp))
                        break;
                    h21 = SB_REAL32_NEGONE;
                    h12 = SB_REAL32_ONE;
                }
            }
        }
    }
    else if f32_ge(tmp, gamsq) {
        if (f32_eq(flag, SB_REAL32_ZERO)) { flag = SB_REAL32_NEGONE; h11 = h22 = SB_REAL32_ONE; }
        else { flag = h21 = SB_REAL32_NEGONE; h12 = SB_REAL32_ONE; }
        if (f32_ge(d2, SB_REAL32_ZERO))
        while (1)
        {
            d2 = f32_mul(d2, rgamsq);
            h21 = f32_mul(h21, gam);
            h22 = f32_mul(h22, gam);
            if (f32_lt(d2, gamsq))
                break;
            h21 = SB_REAL32_NEGONE;
            h12 = SB_REAL32_ONE;
        }
        else /* d2 < SB_REAL32_ZERO */
        {
            tmp = f32_mul(SB_REAL32_NEGONE, gamsq);
            while (1)
            {
                d2 = f32_mul(d2, rgamsq);
                h21 = f32_mul(h21, gam);
                h22 = f32_mul(h22, gam);
                if (f32_gt(d2, tmp))
                    break;
                h21 = SB_REAL32_NEGONE;
                h12 = SB_REAL32_ONE;
            }
        }
    }
   *D1 = d1;
   *D2 = d2;
   *X1 = x1;
   *P = flag;
   if (flag == SB_REAL32_NEGONE) { P[1] = h11; P[2] = h21; P[3] = h12; P[4] = h22; }
   else if (flag == SB_REAL32_ZERO) { P[2] = h21; P[3] = h12; }
   else if (flag == SB_REAL32_ONE) { P[1] = h11; P[4] = h22; }
}
