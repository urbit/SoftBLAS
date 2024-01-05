#include "softblas.h"

uint64_t idamax(uint64_t N, const float64_t *DX, uint64_t incX) {
      float64_t dmax;
      uint64_t i, ix;

      uint64_t idamax = 0;
      if (N < 1 || incX <= 0) return idamax;
      idamax = 1;
      if (N == 1) return idamax;
      if (incX == 1) {
            dmax = f64_abs(DX[0]);
            for (i = 1; i < N; i++) {
                  if (f64_gt(f64_abs(DX[i]), dmax)) {
                        idamax = i;
                        dmax = f64_abs(DX[i]);
                  }
            }
      } else {
            ix = 0;
            dmax = f64_abs(DX[0]);
            ix += incX;
            for (i = 1; i < N; i++) {
                  if (f64_gt(f64_abs(DX[ix]), dmax)) {
                        idamax = i;
                        dmax = f64_abs(DX[ix]);
                  }
                  ix += incX;
            }
      }
      return idamax;
}
