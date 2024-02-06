#include "softblas.h"

float32_t scasum(uint64_t N, const complex32_t *CX, uint64_t incX) {
      float32_t stemp = SB_FLOAT32_ZERO;

      uint64_t NincX = N * incX;
      for (uint64_t i = 0; i < NincX; i += incX) {
         stemp += f32_add(ABS(REAL(CX[i])), ABS(IMAG(CX[i])));
      }

      return stemp;
}
