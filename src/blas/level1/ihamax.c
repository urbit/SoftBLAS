#include "softblas.h"

uint64_t ihamax(uint64_t N, const float16_t *HX, uint64_t incX) {
      float16_t hmax;
      uint64_t i, ix;

      uint64_t ihamax = 0;
      if (N < 1 || incX <= 0) return ihamax;
      ihamax = 1;
      if (N == 1) return ihamax;
      if (incX == 1) {
            hmax = f16_abs(HX[0]);
            for (i = 1; i < N; i++) {
                  if (f16_gt(f16_abs(HX[i]), hmax)) {
                        ihamax = i;
                        hmax = f16_abs(HX[i]);
                  }
            }
      } else {
            ix = 0;
            hmax = f16_abs(HX[0]);
            ix += incX;
            for (i = 1; i < N; i++) {
                  if (f16_gt(f16_abs(HX[ix]), hmax)) {
                        ihamax = i;
                        hmax = f16_abs(HX[ix]);
                  }
                  ix += incX;
            }
      }
      return ihamax;
}
