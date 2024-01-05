#include "softblas.h"

void saxpy(uint64_t N, float32_t SA, float32_t *SX, uint64_t incX, float32_t *SY, uint64_t incY) {
      if (incX == 1 && incY == 1) {
            uint64_t m = N % 4;
            if (m != 0) {
                  for (uint64_t i = 0; i < m; i++) {
                        SY[i] = f32_add(SY[i], f32_mul(SA, SX[i]));
                  }
            }
            if (N < 4) return;
            uint64_t mp1 = m + 1;
            for (uint64_t i = mp1; i < N; i += 4) {
                  SY[i] = f32_add(SY[i], f32_mul(SA, SX[i]));
                  SY[i+1] = f32_add(SY[i+1], f32_mul(SA, SX[i+1]));
                  SY[i+2] = f32_add(SY[i+2], f32_mul(SA, SX[i+2]));
                  SY[i+3] = f32_add(SY[i+3], f32_mul(SA, SX[i+3]));
            }
      } else {
            uint64_t iX = 0;
            uint64_t iY = 0;
            if (incX < 0) iX = (-N + 1) * incX;
            if (incY < 0) iY = (-N + 1) * incY;
            for (uint64_t i = 0; i < N; i++) {
                  SY[iY] = f32_add(SY[iY], f32_mul(SA, SX[iX]));
                  iX += incX;
                  iY += incY;
            }
      }
}
