#include "softblas.h"

void ihrot(const uint64_t N, complex16_t *CX, const int64_t incX, complex16_t *CY, const int64_t incY, const complex16_t c, const complex16_t s, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    complex16_t temp;
    int64_t ix = 0;
    int64_t iy = 0;
    if (incX < 0) ix = (-N + 1) * incX;
    if (incY < 0) iy = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        temp = c16_add(c16_mul(c, CX[ix]), c16_mul(s, CY[iy]));
        CY[iy] = nan_unify_i(c16_sub(c16_mul(c, CY[iy]), c16_mul(s, CX[ix])));
        CX[ix] = nan_unify_i(temp);
        ix += incX;
        iy += incY;
    }
}
