#include "softblas.h"

float32_t snrm2(uint64_t N, const float32_t *X, uint64_t incX) {
      float32_t absXI, norm, scale, ssq;
      
      if (N < 1 || incX < 1) {
            norm = SB_REAL32_ZERO;
      } else if (N == 1) {
            norm = f32_abs(X[0]);
      } else {
            scale = SB_REAL32_ZERO;
            ssq = SB_REAL32_ONE;

            for (uint64_t ix = 0; ix < 1 + (N - 1) * incX; ix += incX) {
                  if (f32_ne(X[ix], SB_REAL32_ZERO)) {
                        absXI = f32_abs(X[ix]);
                        if (f32_lt(scale, absXI)) {
                              ssq = f32_add(SB_REAL32_ONE, f32_mul(ssq, f32_div(f32_mul(scale, scale), f32_mul(absXI, absXI))));
                              scale = absXI;
                        } else {
                              ssq = f32_add(ssq, f32_div(f32_mul(absXI, absXI), f32_mul(scale, scale)));
                        }
                  }
            }

            norm = f32_mul(scale, f32_sqrt(ssq));
      }

      return norm;
}
