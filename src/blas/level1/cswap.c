#include "softblas.h"

//  Swap two complex vectors element-wise: CX <-> CY (BLAS CSWAP).
void cswap(uint64_t N, complex32_t *CX, uint64_t incX, complex32_t *CY, uint64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    complex32_t temp;
    uint64_t iX = 0, iY = 0;
    for (uint64_t i = 0; i < N; i++) {
        temp   = CX[iX];
        CX[iX] = nan_unify_c(CY[iY]);
        CY[iY] = nan_unify_c(temp);
        iX += incX;
        iY += incY;
    }
}
