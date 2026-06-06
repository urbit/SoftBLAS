#include "softblas.h"

complex16_t idotu(uint64_t N, const complex16_t *CX, int64_t incX, const complex16_t *CY, int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    complex16_t idotu = SB_COMPLEX16_ZERO;

    int64_t ix = 0;
    int64_t iy = 0;
    if (incX < 0) ix = (-N + 1) * incX;
    if (incY < 0) iy = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        idotu = c16_add(idotu, c16_mul(CX[ix], CY[iy]));
        ix += incX;
        iy += incY;
    }

    return nan_unify_i(idotu);
}
