#include "softblas.h"

void drotg(
   float64_t *a,    /* INPUT : first rotational elimination parameter */
                    /* OUTPUT: r (see below) */
   float64_t *b,    /* INPUT: second rotational elimination parameter */
                    /* OUTPUT: z (see below) */
   float64_t *c,    /* OUTPUT: cosine */
   float64_t *s     /* OUTPUT: sine */
) {
   float64_t roe, scal, r, z, aa, ab, t0, t1;

   aa = f64_abs(*a);
   ab = f64_abs(*b);
   if f64_gt(aa, ab) roe = *a;
   else roe = *b;
   scal = f64_add(aa, ab);
   if f64_ne(scal, SB_REAL64_ZERO)
   {
      t0 = f64_div(aa, scal); t1 = f64_div(ab, scal);
      r = f64_mul(scal, f64_sqrt(f64_add(f64_mul(t0, t0), f64_mul(t1, t1))));
      if f64_le(roe, SB_REAL64_ZERO) r = f64_neg(r);
      *c = f64_div(*a, r);
      *s = f64_div(*b, r);
      if f64_gt(aa, ab) z = *s;
      else if f64_ne(*c, SB_REAL64_ZERO) z = f64_div(SB_REAL64_ONE, *c);
      else z = SB_REAL64_ONE;
      *a = r;
      *b = z;
   }
   else
   {
      *c = SB_REAL64_ONE;
      *s = *a = *b = SB_REAL64_ZERO;
   }
}
