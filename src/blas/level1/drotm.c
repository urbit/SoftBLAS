#include "softblas.h"

void drotm(const uint64_t N,
           float64_t *X,
           const uint64_t incX,
           float64_t *Y,
           const uint64_t incY,
           const float64_t *P
          ) {
   uint64_t i;
   const float64_t flag = *P;
   float64_t h11, h21, h12, h22, w, z;

   if f64_eq(flag, SB_REAL64_NEGTWO) return;
   if (f64_eq(flag, SB_REAL64_NEGONE))
   {
      h11 = P[1]; h21 = P[2]; h12 = P[3]; h22 = P[4];
      if (incX == 1 && incY == 1)
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X++ = f64_add(f64_mul(w, h11), f64_mul(z, h12));
            *Y++ = f64_add(f64_mul(w, h21), f64_mul(z, h22));
         }
      }
      else
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X = f64_add(f64_mul(w, h11), f64_mul(z, h12));
            *Y = f64_add(f64_mul(w, h21), f64_mul(z, h22));
            X += incX;
            Y += incY;
         }
      }
   }
   else if (f64_eq(flag, SB_REAL64_ZERO))
   {
      h21 = P[2];
      h12 = P[3];
      if (incX == 1 && incY == 1)
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X++ = f64_add(w, f64_mul(z, h12));
            *Y++ = f64_add(f64_mul(w, h21), z);
         }
      }
      else
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X = f64_add(w, f64_mul(z, h12));
            *Y = f64_add(f64_mul(w, h21), z);
            X += incX;
            Y += incY;
         }
      }
   }
   else if (f64_eq(flag, SB_REAL64_ONE))
   {
      h11 = P[1];
      h22 = P[4];
      if (incX == 1 && incY == 1)
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X++ = f64_add(f64_mul(w, h11), z);
            *Y++ = f64_sub(f64_mul(z, h22), w);
         }
      }
      else
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X = f64_add(f64_mul(w, h11), z);
            *Y = f64_sub(f64_mul(z, h22), w);
            X += incX;
            Y += incY;
         }
      }
   }
}
