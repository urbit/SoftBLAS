#include "softblas.h"

void hswap(uint64_t N, float16_t *HX, int64_t incX, float16_t *HY, int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float16_t htemp;

    if (N == 0) return;

    int64_t iX = 0;
    int64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        htemp = HX[iX];
        HX[iX] = nan_unify_h(HY[iY]);
        HY[iY] = nan_unify_h(htemp);
        iX += incX;
        iY += incY;
    }
}
