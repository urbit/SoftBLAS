#include "softblas.h"

void haxpy(uint64_t N, float16_t HA, float16_t *HX, int64_t incX, float16_t *HY, int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    int64_t iX = 0;
    int64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        HY[iY] = f16_add(HY[iY], f16_mul(HA, HX[iX]));
        iX += incX;
        iY += incY;
    }
}
