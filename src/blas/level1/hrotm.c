#include "softblas.h"

void hrotm(const uint16_t N, float16_t *X, const uint16_t incX, float16_t *Y, const uint16_t incY, const float16_t *P) {
   uint16_t i;
   const float16_t flag = *P;
   float16_t h11, h21, h12, h22, w, z;

   if (f16_eq(flag, SB_REAL16_NEGTWO)) return;
   if (f16_eq(flag, SB_REAL16_NEGONE))
   {
      h11 = P[1]; h21 = P[2]; h12 = P[3]; h22 = P[4];
      if (incX == 1 && incY == 1)
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X++ = f16_add(f16_mul(w, h11), f16_mul(z, h12));
            *Y++ = f16_add(f16_mul(w, h21), f16_mul(z, h22));
         }
      }
      else
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X = f16_add(f16_mul(w, h11), f16_mul(z, h12));
            *Y = f16_add(f16_mul(w, h21), f16_mul(z, h22));
            X += incX;
            Y += incY;
         }
      }
   }
   else if (f16_eq(flag, SB_REAL16_ZERO))
   {
      h21 = P[2];
      h12 = P[3];
      if (incX == 1 && incY == 1)
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X++ = f16_add(w, f16_mul(z, h12));
            *Y++ = f16_add(f16_mul(w, h21), z);
         }
      }
      else
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X = f16_add(w, f16_mul(z, h12));
            *Y = f16_add(f16_mul(w, h21), z);
            X += incX;
            Y += incY;
         }
      }
   }
   else if (f16_eq(flag, SB_REAL16_ONE))
   {
      h11 = P[1];
      h22 = P[4];
      if (incX == 1 && incY == 1)
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X++ = f16_add(f16_mul(w, h11), z);
            *Y++ = f16_sub(f16_mul(z, h22), w);
         }
      }
      else
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X = f16_add(f16_mul(w, h11), z);
            *Y = f16_sub(f16_mul(z, h22), w);
            X += incX;
            Y += incY;
         }
      }
   }
}
