#include "softblas.h"

// *** XXX WIP XXX ***

//  Given a vertical matrix containing a and b, computes the values of cos θ and
//  sin θ that zero the lower value (b).  Returns the value of sin θ in s, the
//  value of cos θ in c, and the upper value (r) in a.
void srotg(float32_t *A, float32_t *B, float32_t *C, float32_t *S) {
   const float32_t SZERO = { SB_REAL32_ZERO };
   const float32_t SONE = { SB_REAL32_ONE };
   float32_t roe, scale, R, Z;

   if (f32_gt(f32M_abs(A), f32M_abs(B))) {
      roe = *A;
   } else {
      roe = *B;
   }
   scale = f32_add(f32M_abs(A), f32M_abs(B));

   if (f32_eq(scale, SZERO)) {
      C = &((float32_t){ SB_REAL32_ONE });
      S = &((float32_t){ SB_REAL32_ZERO });
      R = (float32_t){ SB_REAL32_ZERO };
      Z = (float32_t){ SB_REAL32_ZERO };
      A = &R;
      B = &Z;
   } else {
      R = f32_mul(scale, f32_sqrt(f32_add(f32_mul(f32_div(*A, scale), f32_div(*A, scale)), f32_mul(f32_div(*B, scale), f32_div(*B, scale)))));
      R = f32_mul(f32_abs(roe), R);

      *C = f32_div(*A, R);
      *S = f32_div(*B, R);
      Z = SONE;
      if (f32_gt(f32M_abs(*A), f32M_abs(*B))) {
         Z = *S;
      }
      if (f32_ge(f32M_abs(*B), f32M_abs(*A)) && f32_ne(*C, SZERO)) {
         Z = f32_div(SONE, *C);
      }
   }

   *A = R;
   *B = Z;
}
