#include "softblas.h"

float32_t snrm2(uint64_t N, const float32_t *X, uint64_t incX) {
      float32_t ONE = SB_REAL32_ONE;
      float32_t ZERO = SB_REAL32_ZERO;
      float32_t absXI, norm, scale, ssq;
      
      if (N < 1 || incX < 1) {
            norm = ZERO;
      } else if (N == 1) {
            norm = f32_abs(X[0]);
      } else {
            scale = ZERO;
            ssq = ONE;

            for (uint64_t ix = 0; ix < 1 + (N - 1) * incX; ix += incX) {
                  if (!f32_eq(X[ix], ZERO)) {
                        absXI = f32_abs(X[ix]);
                        if (f32_lt(scale, absXI)) {
                              ssq = f32_add(ONE, f32_mul(ssq, f32_div(f32_mul(scale, scale), f32_mul(absXI, absXI))));
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
