#include "softblas.h"

complex16_t idotc(uint64_t N, const complex16_t *CX, int64_t incX, const complex16_t *CY, int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    complex16_t idotc = SB_COMPLEX16_ZERO;

    int64_t ix = 0;
    int64_t iy = 0;
    if (incX < 0) ix = (-N + 1) * incX;
    if (incY < 0) iy = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        idotc = c16_add(idotc, c16_mul(c16_conj(CX[ix]), CY[iy]));
        ix += incX;
        iy += incY;
    }

    return nan_unify_i(idotc);
}
