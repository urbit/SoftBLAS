#include "softblas.h"

void sscal(uint64_t N,
           float16_t HA,
           float16_t *HX,
           uint64_t incX
          ) {
      if (incX == 1) {
            // code for increment equal to 1

            // clean-up loop
            uint64_t m = N % 5;
            if (m != 0) {
                  for (i = 0; i < m; i++) {
                        HX[i] = f16_mul(HA, HX[i]);
                  }
                  if (N < 5) return;
            }
            uint64_t mp1 = m + 1;
            for (uint64_t i = mp1; i < N; i += 5) {
                  HX[i] = f16_mul(HA, HX[i]);
                  HX[i + 1] = f16_mul(HA, HX[i + 1]);
                  HX[i + 2] = f16_mul(HA, HX[i + 2]);
                  HX[i + 3] = f16_mul(HA, HX[i + 3]);
                  HX[i + 4] = f16_mul(HA, HX[i + 4]);
            }
      } else {
            // code for increment not equal to 1
            uint64_t NincX = N * incX;
            for (uint64_t i = 0; i < NincX; i += incX) {
                  HX[i] = f16_mul(HA, HX[i]);
            }
      }
}
