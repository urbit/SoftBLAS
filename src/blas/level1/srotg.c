#include "softblas.h"

void srotg(float32_t *a, float32_t *b, float32_t *c, float32_t *s) {
   float32_t roe, scal, r, z, aa, ab, t0, t1;

   aa = f32_abs(*a);
   ab = f32_abs(*b);
   if f32_gt(aa, ab) roe = *a;
   else roe = *b;
   scal = f32_add(aa, ab);
   if f32_ne(scal, SB_REAL32_ZERO)
   {
      t0 = f32_div(aa, scal); t1 = f32_div(ab, scal);
      r = f32_mul(scal, f32_sqrt(f32_add(f32_mul(t0, t0), f32_mul(t1, t1))));
      if f32_le(roe, SB_REAL32_ZERO) r = f32_neg(r);
      *c = f32_div(*a, r);
      *s = f32_div(*b, r);
      if f32_gt(aa, ab) z = *s;
      else if f32_ne(*c, SB_REAL32_ZERO) z = f32_div(SB_REAL32_ONE, *c);
      else z = SB_REAL32_ONE;
      *a = r;
      *b = z;
   }
   else
   {
      *c = SB_REAL32_ONE;
      *s = *a = *b = SB_REAL32_ZERO;
   }
}
