#include "softblas.h"

void vaxpy(uint64_t N, complex128_t CA, complex128_t *CX, int64_t incX, complex128_t *CY, int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    int64_t iX = 0;
    int64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        CY[iY] = c128_add(CY[iY], c128_mul(CA, CX[iX]));
        nan_unify_v(&(CY[iY]));
        iX += incX;
        iY += incY;
    }
}
