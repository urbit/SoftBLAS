#include "softblas.h"

float32_t sasum(uint64_t N, const float32_t *SX, uint64_t incX) {
      float32_t stemp = SB_REAL32_ZERO;
      
      uint64_t NincX = N * incX;
      for (uint64_t i = 0; i < NincX; i += incX) {
            stemp = f32_add(stemp, f32_abs(SX[i]));
      }

      return stemp;
}
