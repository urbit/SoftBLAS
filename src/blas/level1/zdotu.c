#include "softblas.h"

complex64_t zdotu(uint64_t N, const complex64_t *CX, int64_t incX, const complex64_t *CY, int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    complex64_t zdotu = SB_COMPLEX64_ZERO;

    int64_t ix = 0;
    int64_t iy = 0;
    if (incX < 0) ix = (-N + 1) * incX;
    if (incY < 0) iy = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        zdotu = c64_add(zdotu, c64_mul(CX[ix], CY[iy]));
        ix += incX;
        iy += incY;
    }

    return nan_unify_z(zdotu);
}
