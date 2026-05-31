#include "softblas.h"

void icopy(uint64_t N, const complex16_t *CX, int64_t incX, complex16_t *CY, int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    int64_t ix = 0;
    int64_t iy = 0;
    if (incX < 0) ix = (-N + 1) * incX;
    if (incY < 0) iy = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        CY[iy] = nan_unify_i(CX[ix]);
        ix += incX;
        iy += incY;
    }
}
