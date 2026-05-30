#include "softblas.h"

void sswap(uint64_t N, float32_t *SX, int64_t incX, float32_t *SY, int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float32_t stemp;

    if (N == 0) return;

    int64_t iX = 0;
    int64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        stemp = SX[iX];
        SX[iX] = nan_unify_s(SY[iY]);
        SY[iY] = nan_unify_s(stemp);
        iX += incX;
        iY += incY;
    }
}
