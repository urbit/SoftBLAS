#include "softblas.h"

#include "softblas.h"
#include "softfloat.h"

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
    if (f32_gt(ABS(q1), ABS(q2)))
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
    
    tmp = ABS(d2);
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

void drotmg(float64_t *D1, float64_t *D2, float64_t *X1, const float64_t y1, float64_t *P)
{
    float64_t flag = *P;
    float64_t d1=(*D1), d2=(*D2), x1=(*X1);
    float64_t h11, h21, h12, h22, tmp, u, p1, p2, q1, q2;
    static const float64_t gam=0x40b0000000000000;  // 4096.0
    static const float64_t gamsq = 0x40b0000000000000 * 0x40b0000000000000;
    static const float64_t rgam = SB_REAL64_ONE / 0x40b0000000000000;
    static const float64_t rgamsq = SB_REAL64_ONE /
                                        (0x40b0000000000000 * 0x40b0000000000000);

    if (f64_lt(d1, SB_REAL64_ZERO))
    {
        *P = SB_REAL64_NEGONE;
        *D1 = *D2 = *X1 = P[1] = P[2] = P[3] = P[4] = SB_REAL64_ZERO;
        return;
    }

    p2 = f64_mul(d2, y1);
    if (f64_eq(p2, SB_REAL64_ZERO))
    {
        *P = SB_REAL64_NEGTWO;
        return;
    }

    p1 = f64_mul(d1, x1);
    q2 = f64_mul(p2, y1);
    q1 = f64_mul(p1, x1);
    if (f64_gt(ABS(q1), ABS(q2)))
    {
        h21 = f64_div(f64_mul(SB_REAL64_NEGONE, y1), x1);
        h12 = f64_div(p2, p1);
        u = f64_sub(SB_REAL64_ONE, f64_mul(h12, h21));
        if (f64_gt(u, SB_REAL64_ZERO))
        {
            flag = SB_REAL64_ZERO;
            d1 = f64_div(d1, u);
            d2 = f64_div(d2, u);
            x1 = f64_mul(x1, u);
        }
    }
    else
    {
        if (f64_lt(q2, SB_REAL64_ZERO))
        {
            *P = SB_REAL64_NEGONE;
            *D1 = *D2 = *X1 = P[1] = P[2] = P[3] = P[4] = SB_REAL64_ZERO;
            return;
        }
        flag = SB_REAL64_ONE;
        h11 = f64_div(p1, p2);
        h22 = f64_div(x1, y1);
        u = f64_add(SB_REAL64_ONE, f64_mul(h11, h22));
        tmp = f64_div(d2, u);
        d2 = f64_div(d1, u);
        d1 = tmp;
        x1 = f64_mul(y1, u);

    }

    if (f64_le(d1, rgamsq))
    {
        if (f64_ne(d1, SB_REAL64_ZERO))
        {
            if (f64_eq(flag, SB_REAL64_ZERO)) { flag = SB_REAL64_NEGONE; h11 = h22 = SB_REAL64_ONE; }
            else { flag = h21 = SB_REAL64_NEGONE; h12 = SB_REAL64_ONE; }
            while(1)
            {
                d1 = f64_mul(d1, gamsq);
                x1 = f64_mul(x1, rgam);
                h11 = f64_mul(h11, rgam);
                h12 = f64_mul(h12, rgam);
                if (f64_gt(d1, rgamsq))
                    break;
                h21 = SB_REAL64_NEGONE;
                h12 = SB_REAL64_ONE;
            }
        }
    }
    else if (f64_ge(d1, gamsq))
    {
        if (f64_eq(flag, SB_REAL64_ZERO)) { flag = SB_REAL64_NEGONE; h11 = h22 = SB_REAL64_ONE; }
        else { flag = h21 = SB_REAL64_NEGONE; h12 = SB_REAL64_ONE; }
        while (1)
        {
          d1 = f64_mul(d1, rgamsq);
          x1 = f64_mul(x1, gam);
          h11 = f64_mul(h11, gam);
          h12 = f64_mul(h12, gam);
          if (f64_lt(d1, gamsq))
            break;
          h21 = SB_REAL64_NEGONE;
          h12 = SB_REAL64_ONE;
        }
    }
    
    tmp = ABS(d2);
    if (f64_le(tmp, rgamsq))
    {
        if (!f64_eq(d2, SB_REAL64_ZERO))
        {
            if (f64_eq(flag, SB_REAL64_ZERO)) { flag = SB_REAL64_NEGONE; h11 = h22 = SB_REAL64_ONE; }
            else { flag = h21 = SB_REAL64_NEGONE; h12 = SB_REAL64_ONE; }
            if (f64_ge(d2, SB_REAL64_ZERO))
            while (1)
            {
                d2 = f64_mul(d2, gamsq);
                h21 = f64_mul(h21, rgam);
                h22 = f64_mul(h22, rgam);
                if (f64_gt(d2, rgamsq))
                    break;
                h21 = SB_REAL64_NEGONE;
                h12 = SB_REAL64_ONE;
            }
            else /* d2 < SB_REAL64_ZERO */
            {
                tmp = f64_mul(SB_REAL64_NEGONE, rgamsq);
                while (1)
                {
                    d2 = f64_mul(d2, gamsq);
                    h21 = f64_mul(h21, rgam);
                    h22 = f64_mul(h22, rgam);
                    if (f64_lt(d2, tmp))
                        break;
                    h21 = SB_REAL64_NEGONE;
                    h12 = SB_REAL64_ONE;
                }
            }
        }
    }
    else if f64_ge(tmp, gamsq) {
        if (f64_eq(flag, SB_REAL64_ZERO)) { flag = SB_REAL64_NEGONE; h11 = h22 = SB_REAL64_ONE; }
        else { flag = h21 = SB_REAL64_NEGONE; h12 = SB_REAL64_ONE; }
        if (f64_ge(d2, SB_REAL64_ZERO))
        while (1)
        {
            d2 = f64_mul(d2, rgamsq);
            h21 = f64_mul(h21, gam);
            h22 = f64_mul(h22, gam);
            if (f64_lt(d2, gamsq))
                break;
            h21 = SB_REAL64_NEGONE;
            h12 = SB_REAL64_ONE;
        }
        else /* d2 < SB_REAL64_ZERO */
        {
            tmp = f64_mul(SB_REAL64_NEGONE, gamsq);
            while (1)
            {
                d2 = f64_mul(d2, rgamsq);
                h21 = f64_mul(h21, gam);
                h22 = f64_mul(h22, gam);
                if (f64_gt(d2, tmp))
                    break;
                h21 = SB_REAL64_NEGONE;
                h12 = SB_REAL64_ONE;
            }
        }
    }
   *D1 = d1;
   *D2 = d2;
   *X1 = x1;
   *P = flag;
   if (flag == SB_REAL64_NEGONE) { P[1] = h11; P[2] = h21; P[3] = h12; P[4] = h22; }
   else if (flag == SB_REAL64_ZERO) { P[2] = h21; P[3] = h12; }
   else if (flag == SB_REAL64_ONE) { P[1] = h11; P[4] = h22; }
}

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
    if (f16_gt(ABS(q1), ABS(q2)))
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
    
    tmp = ABS(d2);
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
