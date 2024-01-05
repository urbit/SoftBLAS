#include "softblas.h"

float64_t snrm2(uint64_t N, const float64_t *X, uint64_t incX) {
      float64_t ONE = SB_REAL64_ONE;
      float64_t ZERO = SB_REAL64_ZERO;
      float64_t absXI, norm, scale, ssq;
      
      if (N < 1 || incX < 1) {
            norm = ZERO;
      } else if (N == 1) {
            norm = f64_abs(X[0]);
      } else {
            scale = ZERO;
            ssq = ONE;

            for (uint64_t ix = 0; ix < 1 + (N - 1) * incX; ix += incX) {
                  if (!f64_eq(X[ix], ZERO)) {
                        absXI = f64_abs(X[ix]);
                        if (f64_lt(scale, absXI)) {
                              ssq = f64_add(ONE, f64_mul(ssq, f64_div(f64_mul(scale, scale), f64_mul(absXI, absXI))));
                              scale = absXI;
                        } else {
                              ssq = f64_add(ssq, f64_div(f64_mul(absXI, absXI), f64_mul(scale, scale)));
                        }
                  }
            }

            norm = f64_mul(scale, f64_sqrt(ssq));
      }

      return norm;
}
