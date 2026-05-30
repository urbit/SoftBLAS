#include "softblas.h"

void dswap(uint64_t N, float64_t *DX, int64_t incX, float64_t *DY, int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float64_t dtemp;

    if (N == 0) return;

    int64_t iX = 0;
    int64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        dtemp = DX[iX];
        DX[iX] = nan_unify_d(DY[iY]);
        DY[iY] = nan_unify_d(dtemp);
        iX += incX;
        iY += incY;
    }
}
