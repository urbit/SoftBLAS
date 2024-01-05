#include "softblas.h"

void daxpy(uint64_t N, float64_t HA, float64_t *HX, uint64_t incX, float64_t *HY, uint64_t incY) {
      if (incX == 1 && incY == 1) {
            uint64_t m = N % 4;
            if (m != 0) {
                  for (uint64_t i = 0; i < m; i++) {
                        HY[i] = f64_add(HY[i], f64_mul(HA, HX[i]));
                  }
            }
            if (N < 4) return;
            uint64_t mp1 = m + 1;
            for (uint64_t i = mp1; i < N; i += 4) {
                  HY[i] = f64_add(HY[i], f64_mul(HA, HX[i]));
                  HY[i+1] = f64_add(HY[i+1], f64_mul(HA, HX[i+1]));
                  HY[i+2] = f64_add(HY[i+2], f64_mul(HA, HX[i+2]));
                  HY[i+3] = f64_add(HY[i+3], f64_mul(HA, HX[i+3]));
            }
      } else {
            uint64_t iX = 0;
            uint64_t iY = 0;
            if (incX < 0) iX = (-N + 1) * incX;
            if (incY < 0) iY = (-N + 1) * incY;
            for (uint64_t i = 0; i < N; i++) {
                  HY[iY] = f64_add(HY[iY], f64_mul(HA, HX[iX]));
                  iX += incX;
                  iY += incY;
            }
      }
}
