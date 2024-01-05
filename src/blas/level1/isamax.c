#include "softblas.h"

uint64_t isamax(uint64_t N, const float32_t *SX, uint64_t incX) {
      float32_t smax;
      uint64_t i, ix;

      uint64_t isamax = 0;
      if (N < 1 || incX <= 0) return isamax;
      isamax = 1;
      if (N == 1) return isamax;
      if (incX == 1) {
            smax = f32_abs(SX[0]);
            for (i = 1; i < N; i++) {
                  if (f32_gt(f32_abs(SX[i]), smax)) {
                        isamax = i;
                        smax = f32_abs(SX[i]);
                  }
            }
      } else {
            ix = 0;
            smax = f32_abs(SX[0]);
            ix += incX;
            for (i = 1; i < N; i++) {
                  if (f32_gt(f32_abs(SX[ix]), smax)) {
                        isamax = i;
                        smax = f32_abs(SX[ix]);
                  }
                  ix += incX;
            }
      }
      return isamax;
}
