#include "softblas.h"

void sscal(uint32_t N, float32_t SA, float32_t *SX, uint32_t incX) {
      if (incX == 1) {
            // code for increment equal to 1

            // clean-up loop
            uint64_t m = N % 5;
            if (m != 0) {
                  for (i = 0; i < m; i++) {
                        SX[i] = f32_mul(SA, SX[i]);
                  }
                  if (N < 5) return;
            }
            uint64_t mp1 = m + 1;
            for (uint64_t i = mp1; i < N; i += 5) {
                  SX[i] = f32_mul(SA, SX[i]);
                  SX[i + 1] = f32_mul(SA, SX[i + 1]);
                  SX[i + 2] = f32_mul(SA, SX[i + 2]);
                  SX[i + 3] = f32_mul(SA, SX[i + 3]);
                  SX[i + 4] = f32_mul(SA, SX[i + 4]);
            }
      } else {
            // code for increment not equal to 1
            uint64_t NincX = N * incX;
            for (uint64_t i = 0; i < NincX; i += incX) {
                  SX[i] = f32_mul(SA, SX[i]);
            }
      }
}
