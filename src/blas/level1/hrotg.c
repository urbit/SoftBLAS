#include "softblas.h"

void hrotg(
   float16_t *a,    /* INPUT : first rotational elimination parameter */
                    /* OUTPUT: r (see below) */
   float16_t *b,    /* INPUT: second rotational elimination parameter */
                    /* OUTPUT: z (see below) */
   float16_t *c,    /* OUTPUT: cosine */
   float16_t *s     /* OUTPUT: sine */
) {
   float16_t roe, scal, r, z, aa, ab, t0, t1;

   aa = f16_abs(*a);
   ab = f16_abs(*b);
   if f16_gt(aa, ab) roe = *a;
   else roe = *b;
   scal = f16_add(aa, ab);
   if f16_ne(scal, SB_REAL16_ZERO)
   {
      t0 = f16_div(aa, scal); t1 = f16_div(ab, scal);
      r = f16_mul(scal, f16_sqrt(f16_add(f16_mul(t0, t0), f16_mul(t1, t1))));
      if f16_le(roe, SB_REAL16_ZERO) r = f16_neg(r);
      *c = f16_div(*a, r);
      *s = f16_div(*b, r);
      if f16_gt(aa, ab) z = *s;
      else if f16_ne(*c, SB_REAL16_ZERO) z = f16_div(SB_REAL16_ONE, *c);
      else z = SB_REAL16_ONE;
      *a = r;
      *b = z;
   }
   else
   {
      *c = SB_REAL16_ONE;
      *s = *a = *b = SB_REAL16_ZERO;
   }
}
