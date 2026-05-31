#include "softblas.h"

void iaxpy(uint64_t N, complex16_t CA, complex16_t *CX, int64_t incX, complex16_t *CY, int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
      int64_t iX = 0;
      int64_t iY = 0;
      if (incX < 0) iX = (-N + 1) * incX;
      if (incY < 0) iY = (-N + 1) * incY;
      for (uint64_t i = 0; i < N; i++) {
            CY[iY] = nan_unify_i(c16_add(CY[iY], c16_mul(CA, CX[iX])));
            iX += incX;
            iY += incY;
      }
}
