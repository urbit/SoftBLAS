#include "softblas.h"

//  Swap two complex vectors element-wise: CX <-> CY (BLAS CSWAP).
void zswap(uint64_t N, complex64_t *CX, uint64_t incX, complex64_t *CY, uint64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    complex64_t temp;
    uint64_t iX = 0, iY = 0;
    for (uint64_t i = 0; i < N; i++) {
        temp   = CX[iX];
        CX[iX] = nan_unify_z(CY[iY]);
        CY[iY] = nan_unify_z(temp);
        iX += incX;
        iY += incY;
    }
}
