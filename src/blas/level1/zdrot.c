#include "softblas.h"

void zdrot(const uint64_t N, complex64_t *CX, const int64_t incX, complex64_t *CY, const int64_t incY, const complex64_t c, const complex64_t s, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    complex64_t temp;
    int64_t ix = 0;
    int64_t iy = 0;
    if (incX < 0) ix = (-N + 1) * incX;
    if (incY < 0) iy = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        temp = c64_add(c64_mul(c, CX[ix]), c64_mul(s, CY[iy]));
        CY[iy] = nan_unify_z(c64_sub(c64_mul(c, CY[iy]), c64_mul(s, CX[ix])));
        CX[ix] = nan_unify_z(temp);
        ix += incX;
        iy += incY;
    }
}
