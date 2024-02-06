#include "softblas.h"

float32_t scnrm2(uint64_t N, const complex32_t *CX, uint64_t incX) {
      float32_t norm, scale, ssq, temp;
      
      if (N < 1 || incX < 1) {
            norm = SB_REAL32_ZERO;
      } else {
            scale = SB_REAL32_ZERO;
            ssq = SB_REAL32_ONE;

            for (uint64_t ix = 0; ix < 1 + (N - 1) * incX; ix += incX) {
                  if (f32_ne(REAL(CX[ix]), SB_REAL32_ZERO)) {
                        temp = f32_abs(REAL(CX[ix]));
                        if (f32_lt(scale, temp)) {
                              ssq = f32_add(SB_REAL32_ONE, f32_mul(ssq, f32_div(f32_mul(scale, scale), f32_mul(temp, temp))));
                              scale = temp;
                        } else {
                              ssq = f32_add(ssq, f32_div(f32_mul(temp, temp), f32_mul(scale, scale)));
                        }
                  }
                  if (f32_ne(IMAG(CX[ix]), SB_REAL32_ZERO)) {
                        temp = f32_abs(IMAG(CX[ix]));
                        if (f32_lt(scale, temp)) {
                              ssq = f32_add(SB_REAL32_ONE, f32_mul(ssq, f32_div(f32_mul(scale, scale), f32_mul(temp, temp))));
                              scale = temp;
                        } else {
                              ssq = f32_add(ssq, f32_div(f32_mul(temp, temp), f32_mul(scale, scale)));
                        }
                  }
            }
            norm = f32_mul(scale, f32_sqrt(ssq));
      }

      return norm;
}
