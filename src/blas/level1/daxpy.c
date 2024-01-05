#include "softblas.h"

void daxpy(uint64_t N, float64_t DA, float64_t *DX, uint64_t incX, float64_t *DY, uint64_t incY) {
      if (incX == 1 && incY == 1) {
            uint64_t m = N % 4;
            if (m != 0) {
                  for (uint64_t i = 0; i < m; i++) {
                        DY[i] = f64_add(DY[i], f64_mul(DA, DX[i]));
                  }
            }
            if (N < 4) return;
            uint64_t mp1 = m + 1;
            for (uint64_t i = mp1; i < N; i += 4) {
                  DY[i] = f64_add(DY[i], f64_mul(DA, DX[i]));
                  DY[i+1] = f64_add(DY[i+1], f64_mul(DA, DX[i+1]));
                  DY[i+2] = f64_add(DY[i+2], f64_mul(DA, DX[i+2]));
                  DY[i+3] = f64_add(DY[i+3], f64_mul(DA, DX[i+3]));
            }
      } else {
            uint64_t iX = 0;
            uint64_t iY = 0;
            if (incX < 0) iX = (-N + 1) * incX;
            if (incY < 0) iY = (-N + 1) * incY;
            for (uint64_t i = 0; i < N; i++) {
                  DY[iY] = f64_add(DY[iY], f64_mul(DA, DX[iX]));
                  iX += incX;
                  iY += incY;
            }
      }
}
