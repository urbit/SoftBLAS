#include "softblas.h"

complex32_t cdotu(uint64_t N, const complex32_t *CX, int64_t incX, const complex32_t *CY, int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    complex32_t cdotu = SB_COMPLEX32_ZERO;

    int64_t ix = 0;
    int64_t iy = 0;
    if (incX < 0) ix = (-N + 1) * incX;
    if (incY < 0) iy = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        cdotu = c32_add(cdotu, c32_mul(CX[ix], CY[iy]));
        ix += incX;
        iy += incY;
    }

    return nan_unify_c(cdotu);
}
