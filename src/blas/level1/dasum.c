#include "softblas.h"

float64_t dasum(uint64_t N, const float64_t *DX, uint64_t incX) {
      float64_t dasum = SB_REAL64_ZERO;
      float64_t htemp = SB_REAL64_ZERO;
      
      if (N <= 0 || incX <= 0) {
            return dasum;
      }
      
      if (incX == 1) {
            uint64_t m = N % 6;
            if (m != 0) {
                  for (uint64_t i = 0; i < M; i++) {
                        htemp = f64_add(htemp, f64_abs(DX[i]));
                  }
                  if (N < 6) {
                        dasum = htemp;
                        return dasum;
                  }
            }
            uint64_t mp1 = m + 1;
            for (uint64_t i = mp1; i < N; i += 6) {
                  htemp = f64_add(htemp, f64_add(f64_add(f64_add(f64_add(f64_add(f64_abs(DX[i]), f64_abs(DX[i + 1])), f64_abs(DX[i + 2])), f64_abs(DX[i + 3])), f64_abs(DX[i + 4])), f64_abs(DX[i + 5])));
            }
      } else {
            uint64_t NincX = N * incX;
            for (uint64_t i = 0; i < NincX; i += incX) {
                  htemp = f64_add(htemp, f64_abs(DX[i]));
            }
      }
      dasum = htemp;
      return dasum;
}
