#include "softblas.h"

float16_t hasum(uint64_t N, const float16_t *HX, uint64_t incX) {
      float16_t hasum = SB_REAL16_ZERO;
      float16_t htemp = SB_REAL16_ZERO;
      
      if (N <= 0 || incX <= 0) {
            return hasum;
      }
      
      if (incX == 1) {
            uint64_t m = N % 6;
            if (m != 0) {
                  for (uint64_t i = 0; i < M; i++) {
                        htemp = f16_add(htemp, f16_abs(HX[i]));
                  }
                  if (N < 6) {
                        hasum = htemp;
                        return hasum;
                  }
            }
            uint64_t mp1 = m + 1;
            for (uint64_t i = mp1; i < N; i += 6) {
                  htemp = f16_add(htemp, f16_add(f16_add(f16_add(f16_add(f16_add(f16_abs(HX[i]), f16_abs(HX[i + 1])), f16_abs(HX[i + 2])), f16_abs(HX[i + 3])), f16_abs(HX[i + 4])), f16_abs(HX[i + 5])));
            }
      } else {
            uint64_t NincX = N * incX;
            for (uint64_t i = 0; i < NincX; i += incX) {
                  htemp = f16_add(htemp, f16_abs(HX[i]));
            }
      }
      hasum = htemp;
      return hasum;
}
