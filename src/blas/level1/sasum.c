#include "softblas.h"

float32_t sasum(uint64_t N, const float32_t *SX, uint64_t incX) {
      float32_t sasum = SB_REAL32_ZERO;
      float32_t stemp = SB_REAL32_ZERO;
      
      if (N <= 0 || incX <= 0) {
            return sasum;
      }
      
      if (incX == 1) {
            uint64_t m = N % 6;
            if (m != 0) {
                  for (uint64_t i = 0; i < M; i++) {
                        stemp = f32_add(stemp, f32_abs(SX[i]));
                  }
                  if (N < 6) {
                        sasum = stemp;
                        return sasum;
                  }
            }
            uint64_t mp1 = m + 1;
            for (uint64_t i = mp1; i < N; i += 6) {
                  stemp = f32_add(stemp, f32_add(f32_add(f32_add(f32_add(f32_add(f32_abs(SX[i]), f32_abs(SX[i + 1])), f32_abs(SX[i + 2])), f32_abs(SX[i + 3])), f32_abs(SX[i + 4])), f32_abs(SX[i + 5])));
            }
      } else {
            uint64_t NincX = N * incX;
            for (uint64_t i = 0; i < NincX; i += incX) {
                  stemp = f32_add(stemp, f32_abs(SX[i]));
            }
      }
      sasum = stemp;
      return sasum;
}
