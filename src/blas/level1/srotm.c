#include "softblas.h"

void srotm(const uint32_t N, float32_t *X, const uint32_t incX, float32_t *Y, const uint32_t incY, const float32_t *P) {
   uint64_t i;
   const float32_t flag = *P;
   float32_t h11, h21, h12, h22, w, z;

   if f32_eq(flag, SB_REAL32_NEGTWO) return;
   if (f32_eq(flag, SB_REAL32_NEGONE))
   {
      h11 = P[1]; h21 = P[2]; h12 = P[3]; h22 = P[4];
      if (incX == 1 && incY == 1)
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X++ = f32_add(f32_mul(w, h11), f32_mul(z, h12));
            *Y++ = f32_add(f32_mul(w, h21), f32_mul(z, h22));
         }
      }
      else
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X = f32_add(f32_mul(w, h11), f32_mul(z, h12));
            *Y = f32_add(f32_mul(w, h21), f32_mul(z, h22));
            X += incX;
            Y += incY;
         }
      }
   }
   else if (f32_eq(flag, SB_REAL32_ZERO))
   {
      h21 = P[2];
      h12 = P[3];
      if (incX == 1 && incY == 1)
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X++ = f32_add(w, f32_mul(z, h12));
            *Y++ = f32_add(f32_mul(w, h21), z);
         }
      }
      else
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X = f32_add(w, f32_mul(z, h12));
            *Y = f32_add(f32_mul(w, h21), z);
            X += incX;
            Y += incY;
         }
      }
   }
   else if (f32_eq(flag, SB_REAL32_ONE))
   {
      h11 = P[1];
      h22 = P[4];
      if (incX == 1 && incY == 1)
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X++ = f32_add(f32_mul(w, h11), z);
            *Y++ = f32_sub(f32_mul(z, h22), w);
         }
      }
      else
      {
         for (i = N; i; i--)
         {
            w = *X;
            z = *Y;
            *X = f32_add(f32_mul(w, h11), z);
            *Y = f32_sub(f32_mul(z, h22), w);
            X += incX;
            Y += incY;
         }
      }
   }
}
