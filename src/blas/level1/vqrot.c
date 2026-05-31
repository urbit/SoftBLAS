#include "softblas.h"

//  Apply a rotation [c s; -s c] to a pair of quad-complex vectors (c, s
//  complex, matching csrot).
void vqrot(const uint64_t N, complex128_t *CX, const int64_t incX, complex128_t *CY, const int64_t incY, const complex128_t c, const complex128_t s, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    complex128_t temp;
    int64_t ix = 0;
    int64_t iy = 0;
    if (incX < 0) ix = (-N + 1) * incX;
    if (incY < 0) iy = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        temp = c128_add(c128_mul(c, CX[ix]), c128_mul(s, CY[iy]));
        CY[iy] = c128_sub(c128_mul(c, CY[iy]), c128_mul(s, CX[ix]));
        nan_unify_v(&(CY[iy]));
        CX[ix] = temp;
        nan_unify_v(&(CX[ix]));
        ix += incX;
        iy += incY;
    }
}
