#include "softblas.h"

//  Swap two complex vectors element-wise: CX <-> CY (quad complex).
void vswap(uint64_t N, complex128_t *CX, uint64_t incX, complex128_t *CY, uint64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    complex128_t temp;
    uint64_t iX = 0, iY = 0;
    for (uint64_t i = 0; i < N; i++) {
        temp   = CX[iX];
        CX[iX] = CY[iY];
        nan_unify_v(&(CX[iX]));
        CY[iY] = temp;
        nan_unify_v(&(CY[iY]));
        iX += incX;
        iY += incY;
    }
}
