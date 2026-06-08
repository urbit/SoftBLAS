#include "softblas.h"

complex128_t vdotu(uint64_t N, const complex128_t *CX, int64_t incX, const complex128_t *CY, int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    complex128_t vdotu = SB_COMPLEX128_ZERO;

    int64_t ix = 0;
    int64_t iy = 0;
    if (incX < 0) ix = (-N + 1) * incX;
    if (incY < 0) iy = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        vdotu = c128_add(vdotu, c128_mul(CX[ix], CY[iy]));
        ix += incX;
        iy += incY;
    }

    nan_unify_v(&vdotu);
    return vdotu;
}
