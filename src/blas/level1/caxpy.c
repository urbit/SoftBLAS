#include "softblas.h"

void caxpy(uint64_t N, complex32_t CA, complex32_t *CX, int64_t incX, complex32_t *CY, int64_t incY) {
      int64_t iX = 0;
      int64_t iY = 0;
      if (incX < 0) iX = (-N + 1) * incX;
      if (incY < 0) iY = (-N + 1) * incY;
      for (uint64_t i = 0; i < N; i++) {
            CY[iY] = c32_add(CY[iY], c32_mul(CA, CX[iX]));
            iX += incX;
            iY += incY;
      }
}
