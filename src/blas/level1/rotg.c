#include "softblas.h"
#include "softfloat.h"

/*
 *   | c  s|*|a| = |r|
 *   |-s  c| |b|   |0|
 * This routine returns:
 *    r = sigma * sqrt(a^2 + b^2), where
 *      sigma = sign(a) if abs(a) > abs(b)
 *      sigma = sign(b) if abs(a) <= abs(b)
 *    r is returned in *a
 *
 *    z = s     if (abs(a) > abs(b))
 *    z = 1/c   if (abs(a) <= abs(b) && c != 0 && r != 0)
 *    z = 1     if (abs(a) <= abs(b) && c == 0 && r != 0)
 *    z = 0     if (r == 0)
 *    z is returned in *b
 *
 *    c : cosign of the angle of (Givens) rotation
 *    c = a/r   if (r != 0)
 *    c = 1     if (r == 0)
 *
 *    s : sine of the angle of (Givens) rotation
 *    s = b/r   if (r != 0)
 *    s = 0     if (r == 0)
 *    FURTHER DETAILS:
 *       http://publib.boulder.ibm.com/infocenter/clresctr/vxrx/index.jsp?topic=/com.ibm.cluster.essl43.guideref.doc/am501_hsrotg.html
 *
 */

void srotg(
   float32_t *a,    /* INPUT : first rotational elimination parameter */
                    /* OUTPUT: r (see below) */
   float32_t *b,    /* INPUT: second rotational elimination parameter */
                    /* OUTPUT: z (see below) */
   float32_t *c,    /* OUTPUT: cosine */
   float32_t *s     /* OUTPUT: sine */
) {
   float32_t roe, scal, r, z, aa, ab, t0, t1;

   aa = ABS(*a);
   ab = ABS(*b);
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

void drotg(
   float64_t *a,    /* INPUT : first rotational elimination parameter */
                    /* OUTPUT: r (see below) */
   float64_t *b,    /* INPUT: second rotational elimination parameter */
                    /* OUTPUT: z (see below) */
   float64_t *c,    /* OUTPUT: cosine */
   float64_t *s     /* OUTPUT: sine */
) {
   float64_t roe, scal, r, z, aa, ab, t0, t1;

   aa = ABS(*a);
   ab = ABS(*b);
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

#include "softblas.h"
#include "softfloat.h"

void hrotg(
   float16_t *a,    /* INPUT : first rotational elimination parameter */
                    /* OUTPUT: r (see below) */
   float16_t *b,    /* INPUT: second rotational elimination parameter */
                    /* OUTPUT: z (see below) */
   float16_t *c,    /* OUTPUT: cosine */
   float16_t *s     /* OUTPUT: sine */
) {
   float16_t roe, scal, r, z, aa, ab, t0, t1;

   aa = ABS(*a);
   ab = ABS(*b);
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
