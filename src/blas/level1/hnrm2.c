#include "softblas.h"

float16_t hnrm2(uint64_t N, const float16_t *X, uint64_t incX) {
      float16_t ONE = SB_REAL16_ONE;
      float16_t ZERO = SB_REAL16_ZERO;
      float16_t absXI, norm, scale, ssq;
      
      if (N < 1 || incX < 1) {
            norm = ZERO;
      } else if (N == 1) {
            norm = f16_abs(X[0]);
      } else {
            scale = ZERO;
            ssq = ONE;

            for (uint64_t ix = 0; ix < 1 + (N - 1) * incX; ix += incX) {
                  if (!f16_eq(X[ix], ZERO)) {
                        absXI = f16_abs(X[ix]);
                        if (f16_lt(scale, absXI)) {
                              ssq = f16_add(ONE, f16_mul(ssq, f16_div(f16_mul(scale, scale), f16_mul(absXI, absXI))));
                              scale = absXI;
                        } else {
                              ssq = f16_add(ssq, f16_div(f16_mul(absXI, absXI), f16_mul(scale, scale)));
                        }
                  }
            }

            norm = f16_mul(scale, f16_sqrt(ssq));
      }

      return norm;
}
